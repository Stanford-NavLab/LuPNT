#include <lupnt/conversions/frame_conversions.h>
#include <lupnt/conversions/frame_converter.h>
#include <lupnt/conversions/time_conversions.h>
#include <lupnt/interfaces/spice.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // Assert R is orthonormal with determinant +1 (a proper rotation).
  void RequireProperRotation(const Mat3& R, double tol) {
    Mat3 I = R * R.transpose();
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) REQUIRE_THAT(I(i, j).val(), WithinAbs(i == j ? 1.0 : 0.0, tol));
    REQUIRE_THAT(R.determinant().val(), WithinAbs(1.0, tol));
  }
  void RequireProperRotation(const Mat3d& R, double tol) {
    Mat3d I = R * R.transpose();
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) REQUIRE_THAT(I(i, j), WithinAbs(i == j ? 1.0 : 0.0, tol));
    REQUIRE_THAT(R.determinant(), WithinAbs(1.0, tol));
  }

  const Vec6 kRvSample = (Vec6() << 1.7e6, -2.3e6, 3.1e6, 850.0, -1200.0, 640.0).finished();
}  // namespace

// ---------------------------------------------------------------------------
// MoonPa <-> MoonMe : a constant IAU frame-bias rotation (no epoch, no SPICE).
// ---------------------------------------------------------------------------
TEST_CASE("conversions.frame_conversions_extra.moon_pa_me") {
  SECTION("RotMoonPaToMe is a proper rotation") { RequireProperRotation(RotMoonPaToMe(), 1e-14); }

  SECTION("MoonPaToMe / MoonMeToPa round-trip position and velocity") {
    Vec6 rv_pa = kRvSample;
    Vec6 rt = MoonMeToPa(MoonPaToMe(rv_pa));
    for (int i = 0; i < 6; ++i) REQUIRE_THAT(rt(i).val(), WithinAbs(rv_pa(i).val(), 1e-8));
  }

  SECTION("Vec3 and Vec6 overloads agree on the position block") {
    Vec6 rv_pa = kRvSample;
    Vec3 r_me_v3 = MoonPaToMe(Vec3(rv_pa.head(3)));
    Vec6 rv_me_v6 = MoonPaToMe(rv_pa);
    for (int i = 0; i < 3; ++i) REQUIRE_THAT(rv_me_v6(i).val(), WithinAbs(r_me_v3(i).val(), 1e-12));

    Vec3 r_pa_v3 = MoonMeToPa(Vec3(rv_me_v6.head(3)));
    for (int i = 0; i < 3; ++i) REQUIRE_THAT(r_pa_v3(i).val(), WithinAbs(rv_pa(i).val(), 1e-8));
  }

  SECTION("the bias rotation preserves vector norms") {
    Vec3 r_pa(1.2e6, -3.4e5, 9.8e5);
    Vec3 r_me = MoonPaToMe(r_pa);
    REQUIRE_THAT(r_me.norm().val(), WithinRel(r_pa.norm().val(), 1e-12));
  }
}

// ---------------------------------------------------------------------------
// GCRF <-> EME2000 constant frame bias: exact vs. first/second-order forms.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.frame_conversions_extra.gcrf_eme_orders") {
  SECTION("exact GCRF->EME bias is a proper rotation") {
    RequireProperRotation(RotGcrfToEme(), 1e-13);
  }

  SECTION("first- and second-order forms converge toward the exact rotation") {
    Mat3d R_exact = RotGcrfToEme();
    Mat3d R1 = RotGcrfToEmeFirstOrder();
    Mat3d R2 = RotGcrfToEmeSecondOrder();

    double err1 = (R1 - R_exact).cwiseAbs().maxCoeff();
    double err2 = (R2 - R_exact).cwiseAbs().maxCoeff();

    // The bias angles are ~tens of milliarcsec (~1e-7 rad), so both truncations
    // are extremely close to exact, and the second order is at least as good.
    REQUIRE(err1 < 1e-13);
    REQUIRE(err2 <= err1);
  }

  SECTION("Vec3 and Vec6 GcrfToEme overloads agree on the position block") {
    Vec6 rv_gcrf = kRvSample;
    Vec3 r_v3 = GcrfToEme(Vec3(rv_gcrf.head(3)));
    Vec6 rv_v6 = GcrfToEme(rv_gcrf);
    for (int i = 0; i < 3; ++i) REQUIRE_THAT(rv_v6(i).val(), WithinAbs(r_v3(i).val(), 1e-9));
  }
}

// ---------------------------------------------------------------------------
// Translational hub conversions (bundled DE ephemeris; no SPICE kernels
// required). These are pure additions/subtractions of a body ephemeris state.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.frame_conversions_extra.translations") {
  std::vector<double> epochs = {0.0, 5.0 * 86400.0, -7.0 * 86400.0};

  SECTION("GcrfToIcrf / IcrfToGcrf round-trip and Vec3/Vec6 agree") {
    for (double t : epochs) {
      INFO("t_tdb = " << t);
      Real t_tdb(t);
      Vec6 rv = kRvSample;
      Vec6 rt = IcrfToGcrf(t_tdb, GcrfToIcrf(t_tdb, rv));
      for (int i = 0; i < 6; ++i) REQUIRE_THAT(rt(i).val(), WithinAbs(rv(i).val(), 1e-6));

      Vec3 r_v3 = GcrfToIcrf(t_tdb, Vec3(rv.head(3)));
      Vec6 r_v6 = GcrfToIcrf(t_tdb, rv);
      for (int i = 0; i < 3; ++i) REQUIRE_THAT(r_v6(i).val(), WithinAbs(r_v3(i).val(), 1e-6));
    }
  }

  SECTION("GcrfToMoonCi / MoonCiToGcrf round-trip; zero input recovers Earth->Moon vector") {
    for (double t : epochs) {
      INFO("t_tdb = " << t);
      Real t_tdb(t);
      Vec6 rv = kRvSample;
      Vec6 rt = MoonCiToGcrf(t_tdb, GcrfToMoonCi(t_tdb, rv));
      for (int i = 0; i < 6; ++i) REQUIRE_THAT(rt(i).val(), WithinAbs(rv(i).val(), 1e-6));

      // GcrfToMoonCi(0) = -(Earth->Moon), whose magnitude is the lunar distance.
      Vec3 offset = GcrfToMoonCi(t_tdb, Vec3(0.0, 0.0, 0.0));
      REQUIRE(offset.norm().val() > 3.5e8);  // perigee ~3.63e8 m
      REQUIRE(offset.norm().val() < 4.1e8);  // apogee ~4.05e8 m
    }
  }
}

// ---------------------------------------------------------------------------
// MoonCi <-> MoonPa analytic path (DE Chebyshev lunar mantle libration; bundled
// data, no InitFrameConversionFromSpice fit installed).
// ---------------------------------------------------------------------------
TEST_CASE("conversions.frame_conversions_extra.moon_ci_pa_analytic") {
  ClearFrameConversionFit();  // ensure the analytic (GetLunarMantleData) path
  std::vector<double> epochs = {0.0, 2.0 * 86400.0, 5.0 * 86400.0};

  SECTION("RotMoonCiToPa is a proper rotation") {
    for (double t : epochs) {
      INFO("t_tdb = " << t);
      RequireProperRotation(RotMoonCiToPa(Real(t)), 1e-12);
    }
  }

  SECTION("MoonCiToPa / MoonPaToCi round-trip position and velocity") {
    for (double t : epochs) {
      INFO("t_tdb = " << t);
      Real t_tdb(t);
      Vec6 rv = kRvSample;
      Vec6 rt = MoonPaToCi(t_tdb, MoonCiToPa(t_tdb, rv));
      for (int i = 0; i < 6; ++i) REQUIRE_THAT(rt(i).val(), WithinAbs(rv(i).val(), 1e-6));

      Vec3 r_v3 = MoonCiToPa(t_tdb, Vec3(rv.head(3)));
      Vec6 r_v6 = MoonCiToPa(t_tdb, rv);
      for (int i = 0; i < 3; ++i) REQUIRE_THAT(r_v6(i).val(), WithinAbs(r_v3(i).val(), 1e-9));
    }
  }

  SECTION("analytic R_dot matches a central finite difference of RotMoonCiToPa") {
    Real t_tdb(2.0 * 86400.0);
    double dt = 1.0;
    Mat3 R_dot;
    RotMoonCiToPa(t_tdb, &R_dot);
    Mat3 R_plus = RotMoonCiToPa(t_tdb + dt);
    Mat3 R_minus = RotMoonCiToPa(t_tdb - dt);
    Mat3 R_dot_num = (R_plus - R_minus) / (2.0 * dt);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        REQUIRE_THAT(R_dot(i, j).val(), WithinAbs(R_dot_num(i, j).val(), 1e-11));
  }
}

// ---------------------------------------------------------------------------
// MoonCi <-> MoonOp (Ely frozen-orbit "orbit plane" frame; bundled ephemeris).
// ---------------------------------------------------------------------------
TEST_CASE("conversions.frame_conversions_extra.moon_op") {
  std::vector<double> epochs = {0.0, 5.0 * 86400.0};

  SECTION("RotMoonOpToCi is a proper rotation") {
    for (double t : epochs) {
      INFO("t_tdb = " << t);
      RequireProperRotation(RotMoonOpToCi(Real(t)), 1e-12);
    }
  }

  SECTION("MoonCiToOp / MoonOpToCi round-trip position and velocity") {
    for (double t : epochs) {
      INFO("t_tdb = " << t);
      Real t_tdb(t);
      Vec6 rv = kRvSample;
      Vec6 rt = MoonOpToCi(t_tdb, MoonCiToOp(t_tdb, rv));
      for (int i = 0; i < 6; ++i) REQUIRE_THAT(rt(i).val(), WithinAbs(rv(i).val(), 1e-8));

      Vec3 r_v3 = MoonCiToOp(t_tdb, Vec3(rv.head(3)));
      Vec6 r_v6 = MoonCiToOp(t_tdb, rv);
      for (int i = 0; i < 3; ++i) REQUIRE_THAT(r_v6(i).val(), WithinAbs(r_v3(i).val(), 1e-9));
    }
  }
}

// ---------------------------------------------------------------------------
// GCRF <-> EMR (Earth-Moon rotating/synodic frame; bundled ephemeris).
// ---------------------------------------------------------------------------
TEST_CASE("conversions.frame_conversions_extra.emr") {
  std::vector<double> epochs = {0.0, 5.0 * 86400.0};

  SECTION("GcrfToEmr / EmrToGcrf round-trip full state") {
    for (double t : epochs) {
      INFO("t_tdb = " << t);
      Real t_tdb(t);
      Vec6 rv = kRvSample;
      Vec6 rt = EmrToGcrf(t_tdb, GcrfToEmr(t_tdb, rv));
      for (int i = 0; i < 6; ++i) REQUIRE_THAT(rt(i).val(), WithinAbs(rv(i).val(), 1e-6));
    }
  }

  SECTION("Vec3 (position-only) round-trip recovers the input position") {
    for (double t : epochs) {
      INFO("t_tdb = " << t);
      Real t_tdb(t);
      Vec3 r = kRvSample.head(3);
      Vec3 rt = EmrToGcrf(t_tdb, GcrfToEmr(t_tdb, r));
      for (int i = 0; i < 3; ++i) REQUIRE_THAT(rt(i).val(), WithinAbs(r(i).val(), 1e-6));
    }
  }
}

// ---------------------------------------------------------------------------
// Generic-IAU planet body orientation (Mercury..Neptune; bundled ephemeris).
// ---------------------------------------------------------------------------
TEST_CASE("conversions.frame_conversions_extra.planet_iau") {
  SECTION("HasIauOrientation covers the planets but not Earth/Moon/Sun") {
    REQUIRE(HasIauOrientation(BodyId::MERCURY));
    REQUIRE(HasIauOrientation(BodyId::MARS));
    REQUIRE(HasIauOrientation(BodyId::JUPITER));
    REQUIRE(HasIauOrientation(BodyId::NEPTUNE));
    REQUIRE_FALSE(HasIauOrientation(BodyId::EARTH));
    REQUIRE_FALSE(HasIauOrientation(BodyId::MOON));
    REQUIRE_FALSE(HasIauOrientation(BodyId::SUN));
  }

  SECTION("RotBodyCiToFixed is a proper rotation and its R_dot matches a finite difference") {
    Real t_tdb(5.0 * 86400.0);
    for (BodyId body : {BodyId::MARS, BodyId::JUPITER}) {
      INFO("body = " << static_cast<int>(body));
      Mat3 R_dot;
      Mat3 R = RotBodyCiToFixed(t_tdb, body, &R_dot);
      RequireProperRotation(R, 1e-12);

      double dt = 1.0;
      Mat3 R_plus = RotBodyCiToFixed(t_tdb + dt, body);
      Mat3 R_minus = RotBodyCiToFixed(t_tdb - dt, body);
      Mat3 R_dot_num = (R_plus - R_minus) / (2.0 * dt);
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          REQUIRE_THAT(R_dot(i, j).val(), WithinAbs(R_dot_num(i, j).val(), 1e-9));
    }
  }

  SECTION("BodyCiToFixed / BodyFixedToCi round-trip position and velocity") {
    Real t_tdb(5.0 * 86400.0);
    Vec6 rv = kRvSample;
    Vec6 rt = BodyFixedToCi(t_tdb, BodyCiToFixed(t_tdb, rv, BodyId::MARS), BodyId::MARS);
    for (int i = 0; i < 6; ++i) REQUIRE_THAT(rt(i).val(), WithinAbs(rv(i).val(), 1e-6));

    Vec3 r_v3 = BodyCiToFixed(t_tdb, Vec3(rv.head(3)), BodyId::MARS);
    Vec6 r_v6 = BodyCiToFixed(t_tdb, rv, BodyId::MARS);
    for (int i = 0; i < 3; ++i) REQUIRE_THAT(r_v6(i).val(), WithinAbs(r_v3(i).val(), 1e-9));
  }

  SECTION(
      "PlanetCiToIcrf / IcrfToPlanetCi round-trip; zero input recovers the SSB->planet vector") {
    Real t_tdb(0.0);
    Vec6 rv = kRvSample;
    Vec6 rt = IcrfToPlanetCi(t_tdb, PlanetCiToIcrf(t_tdb, rv, BodyId::MARS), BodyId::MARS);
    for (int i = 0; i < 6; ++i) REQUIRE_THAT(rt(i).val(), WithinAbs(rv(i).val(), 1e-4));

    // PlanetCiToIcrf(0) = SSB->planet, of order the planet's heliocentric range.
    Vec3 offset = PlanetCiToIcrf(t_tdb, Vec3(0.0, 0.0, 0.0), BodyId::MARS);
    REQUIRE(offset.norm().val() > 1.5e11);  // Mars is >1 AU from the SSB
  }
}

// The TIO locator s' (IERS Conventions 2010, Eq. 5.13) is s' = -47 uas * T
// with T in JULIAN CENTURIES of TT. LuPNT holds times as seconds from J2000,
// so the conversion must divide by DAYS_CENTURY * SECS_DAY. Dividing by
// DAYS_CENTURY alone inflated s' by 86400x -- ~5e-6 rad by the mid-2020s,
// which is ~31 m of spurious rotation about the polar axis at the equator.
//
// s' is not exposed directly, so it is recovered from the polar-motion matrix:
// R_po = RotX(-yp) * RotY(-xp) * RotZ(sp), whose (0,1) entry is sp to first
// order in the three small angles.
TEST_CASE("conversions.frame_conversions_extra.tio_locator_magnitude") {
  spice::LoadSpiceKernel();

  // ~2025: T ~ 0.25 Julian centuries, so |s'| ~ 47e-6 * 0.25 arcsec ~ 6e-11 rad.
  Real t_tdb = GregorianToTime(2025, 1, 1, 0, 0, 0.0);
  Mat3 R_po = RotPolarMotion(t_tdb);

  double sp = R_po(0, 1).val();
  INFO("recovered s' = " << sp << " rad");

  // Correct magnitude is ~6e-11 rad; the defect gave ~5e-6 rad. A 1e-8 bound
  // separates the two by orders of magnitude without over-fitting the value.
  REQUIRE(std::abs(sp) < 1e-8);

  // s' grows linearly with T and must stay negative (the -47 uas coefficient).
  Real t_2075 = GregorianToTime(2075, 1, 1, 0, 0, 0.0);
  double sp_2075 = RotPolarMotion(t_2075)(0, 1).val();
  INFO("s'(2075) = " << sp_2075 << " rad");
  REQUIRE(sp_2075 < sp);  // more negative later
  REQUIRE(std::abs(sp_2075) < 1e-8);
}
