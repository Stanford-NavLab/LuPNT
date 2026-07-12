#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "lupnt/environment/plasma/env/orbit_utils.h"
#include "lupnt/environment/plasma/env/time_utils.h"

using namespace pecsim;
using namespace Catch::Matchers;

// ------------------------------- time_utils --------------------------------

TEST_CASE("environment.plasma.time_utils.doy_mmdd_roundtrip") {
  SECTION("non-leap year day 60 is March 1") {
    int mm, dd;
    doy_to_mmdd(2023, 60, mm, dd);
    REQUIRE(mm == 3);
    REQUIRE(dd == 1);
  }
  SECTION("leap year day 60 is February 29") {
    int mm, dd;
    doy_to_mmdd(2024, 60, mm, dd);
    REQUIRE(mm == 2);
    REQUIRE(dd == 29);
  }
  SECTION("mmdd_to_doy inverts doy_to_mmdd") {
    for (int doy = 1; doy <= 365; ++doy) {
      int mm, dd, doy_back;
      doy_to_mmdd(2023, doy, mm, dd);
      mmdd_to_doy(2023, mm, dd, doy_back);
      REQUIRE(doy_back == doy);
    }
  }
  SECTION("out-of-range day throws") {
    int mm, dd;
    REQUIRE_THROWS(doy_to_mmdd(2023, 0, mm, dd));
    REQUIRE_THROWS(doy_to_mmdd(2023, 400, mm, dd));
  }
}

TEST_CASE("environment.plasma.time_utils.itime_roundtrip") {
  DateTime dt{2024, 100, 13, 45, 30.0};
  std::array<int, 2> itime = datetime_to_itime(dt);
  REQUIRE(itime[0] == 2024 * 1000 + 100);
  REQUIRE(itime[1] == (13 * 3600 + 45 * 60 + 30) * 1000);

  DateTime back = itime_to_datetime(itime);
  REQUIRE(back.year == 2024);
  REQUIRE(back.doy == 100);
  REQUIRE(back.hour == 13);
  REQUIRE(back.min == 45);
  REQUIRE_THAT(back.sec, WithinAbs(30.0, 1.0e-9));
}

TEST_CASE("environment.plasma.time_utils.mjd_conversions") {
  SECTION("gregorian_to_mjd of J2000 midnight is 51544") {
    REQUIRE_THAT(gregorian_to_mjd(2000, 1, 1, 0, 0, 0.0), WithinAbs(51544.0, 1.0e-9));
  }
  SECTION("datetime_to_mjd of 2000 doy 1 is 51544") {
    DateTime dt{2000, 1, 0, 0, 0.0};
    REQUIRE_THAT(datetime_to_mjd(dt), WithinAbs(51544.0, 1.0e-6));
  }
  SECTION("mjd_to_datetime inverts datetime_to_mjd") {
    DateTime dt{2021, 200, 6, 30, 0.0};
    double mjd = datetime_to_mjd(dt);
    DateTime back = mjd_to_datetime(mjd);
    REQUIRE(back.year == 2021);
    REQUIRE(back.doy == 200);
    REQUIRE(back.hour == 6);
    REQUIRE(back.min == 30);
  }
  SECTION("tj2000 <-> mjd is an exact inverse") {
    double mjd = 59000.25;
    double t = mjd_to_tj2000(mjd);
    REQUIRE_THAT(tj2000_to_mjd(t), WithinAbs(mjd, 1.0e-9));
  }
}

TEST_CASE("environment.plasma.time_utils.local_time_roundtrip") {
  for (double mlt : {0.0, 6.0, 12.0, 18.0, 23.5}) {
    double along = lt_to_long(mlt);
    double mlt_back = long_to_lt(along);
    REQUIRE_THAT(mlt_back, WithinAbs(mlt, 1.0e-9));
  }
}

// ------------------------------- orbit_utils -------------------------------

TEST_CASE("environment.plasma.orbit_utils.rotation_matrices") {
  SECTION("RotXd(0) is the identity") {
    Mat3d R = RotXd(0.0);
    REQUIRE_THAT((R - Mat3d::Identity()).norm(), WithinAbs(0.0, 1.0e-14));
  }
  SECTION("RotXd is orthonormal with unit determinant") {
    double a = 0.7;
    Mat3d R = RotXd(a);
    REQUIRE_THAT((R.transpose() * R - Mat3d::Identity()).norm(), WithinAbs(0.0, 1.0e-13));
    REQUIRE_THAT(R.determinant(), WithinAbs(1.0, 1.0e-13));
    // Rotation about x leaves the x-axis fixed.
    Vec3d ex(1.0, 0.0, 0.0);
    REQUIRE_THAT((R * ex - ex).norm(), WithinAbs(0.0, 1.0e-14));
  }
}

TEST_CASE("environment.plasma.orbit_utils.anomaly_conversions") {
  const double e = 0.3;
  SECTION("ecc <-> mean round trip") {
    for (double E : {0.1, 1.0, 2.5, 4.0, 6.0}) {
      double M = ecc2mean(E, e);
      double E_back = mean2ecc(M, e);
      REQUIRE_THAT(E_back, WithinAbs(E, 1.0e-8));
    }
  }
  SECTION("true <-> ecc round trip") {
    for (double nu : {0.2, 1.2, 2.0}) {
      double E = true2ecc(nu, e);
      double nu_back = ecc2true(E, e);
      REQUIRE_THAT(nu_back, WithinAbs(nu, 1.0e-9));
    }
  }
  SECTION("mean2true is consistent with mean2ecc + ecc2true") {
    double M = 1.1;
    double E = mean2ecc(M, e);
    double nu_direct = mean2true(M, e);
    double nu_composed = ecc2true(E, e);
    // Both are wrapped to [0, 2pi); compare via sin/cos to avoid branch-cut issues.
    REQUIRE_THAT(std::sin(nu_direct), WithinAbs(std::sin(nu_composed), 1.0e-9));
    REQUIRE_THAT(std::cos(nu_direct), WithinAbs(std::cos(nu_composed), 1.0e-9));
  }
}

TEST_CASE("environment.plasma.orbit_utils.coe_cart_roundtrip") {
  // Earth GM in km^3/s^2 (the module works in km); a in km.
  const double GM = 3.986004418e5;
  Vec6d coe;
  coe << 7000.0, 0.05, 0.9, 1.2, 0.5, 0.8;  // a[km], e, i, Omega, omega, M [rad]

  Vec6d rv = coe2cart(coe, GM);
  // Sanity: perigee/apogee bound the radius.
  double r = rv.head<3>().norm();
  REQUIRE(r > coe[0] * (1.0 - coe[1]) - 1.0);
  REQUIRE(r < coe[0] * (1.0 + coe[1]) + 1.0);

  Vec6d coe_back = cart2coe(rv, GM);
  REQUIRE_THAT(coe_back[0], WithinAbs(coe[0], 1.0e-6));  // a
  REQUIRE_THAT(coe_back[1], WithinAbs(coe[1], 1.0e-9));  // e
  REQUIRE_THAT(coe_back[2], WithinAbs(coe[2], 1.0e-9));  // i
  REQUIRE_THAT(coe_back[3], WithinAbs(coe[3], 1.0e-9));  // Omega
  REQUIRE_THAT(coe_back[4], WithinAbs(coe[4], 1.0e-9));  // omega
  REQUIRE_THAT(coe_back[5], WithinAbs(coe[5], 1.0e-9));  // M
}

TEST_CASE("environment.plasma.orbit_utils.propagate_coe") {
  const double GM = 3.986004418e5;
  Vec6d coe;
  coe << 7000.0, 0.01, 0.5, 0.0, 0.0, 0.0;
  double n = std::sqrt(GM / std::pow(coe[0], 3));  // mean motion [rad/s]
  double period = 2.0 * M_PI / n;

  // A full period returns the mean anomaly (and thus all elements) to the start.
  Vec6d after = propagate_coe(coe, GM, period);
  REQUIRE_THAT(std::sin(after[5]), WithinAbs(std::sin(coe[5]), 1.0e-6));
  REQUIRE_THAT(std::cos(after[5]), WithinAbs(std::cos(coe[5]), 1.0e-6));
  // Non-anomaly elements are unchanged.
  REQUIRE_THAT(after[0], WithinAbs(coe[0], 1.0e-9));
  REQUIRE_THAT(after[2], WithinAbs(coe[2], 1.0e-12));
}
