#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <filesystem>
#include <vector>

#include "lupnt/core/file.h"
#include "lupnt/environment/plasma/env/time_utils.h"
#include "lupnt/environment/plasma/gcpm/iri_interface.h"
#include "lupnt/environment/plasma/tec/raytrace.h"

using namespace pecsim;
using namespace Catch::Matchers;

// Exercises the plasma ray-tracing electron-density entry point compute_ne and the
// pure-geometry helpers azel_to_unitvec / unitvec_to_azel / is_method_rk4.
//
// Robustness: compute_ne is driven through the *Fortran* GCPM (use_fortran_gcpm =
// true), the physically-validated path. The pure-C++ GCPM has a known NaN defect on
// the low-latitude branch (see test_plasma_gcpm.cc), so we deliberately avoid it and
// assert finiteness / non-negativity / radial monotonicity rather than exact values.

namespace {
  struct CwdGuard {
    std::filesystem::path saved;
    CwdGuard() : saved(std::filesystem::current_path()) {}
    ~CwdGuard() { std::filesystem::current_path(saved); }
  };

  void SetupPlasma() {
    std::filesystem::path plasma_base = lupnt::GetDataPath() / "plasma";
    set_base_path(plasma_base.string());
    set_iri_model(IRIModel::IRI_2020);
  }

  // A daytime, mid-solar-cycle epoch inside the bundled apf107 / ig_rz range,
  // expressed as seconds since J2000 (what compute_ne expects).
  double SampleTj2000() {
    DateTime dt{2015, 172, 12, 0, 0.0};
    return mjd_to_tj2000(datetime_to_mjd(dt));
  }

  RayTraceConfig FortranConfig() {
    RayTraceConfig config;
    config.freq_Hz = freq_L1;
    config.kp = 3.0;
    config.use_fortran_gcpm = true;
    return config;
  }
}  // namespace

TEST_CASE("environment.plasma.raytrace.geometry_helpers") {
  SECTION("is_method_rk4 recognizes RK4 spellings only") {
    REQUIRE(is_method_rk4("RK4"));
    REQUIRE(is_method_rk4("rk4"));
    REQUIRE_FALSE(is_method_rk4("Euler"));
    REQUIRE_FALSE(is_method_rk4("euler"));
    REQUIRE_FALSE(is_method_rk4(""));
  }

  SECTION("azel_to_unitvec returns a unit vector with the right octant") {
    Vec3d up = azel_to_unitvec(0.0, M_PI / 2.0);  // straight up (+z)
    REQUIRE_THAT(up.norm(), WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(up[2], WithinAbs(1.0, 1e-12));

    Vec3d east = azel_to_unitvec(M_PI / 2.0, 0.0);  // az=90 deg, horizon (+y)
    REQUIRE_THAT(east.norm(), WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(east[1], WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(east[0], WithinAbs(0.0, 1e-12));
  }

  SECTION("azel<->unitvec round-trips") {
    const double azs[] = {-2.5, -1.0, 0.0, 0.7, 2.9};
    const double els[] = {-1.2, -0.3, 0.0, 0.5, 1.3};
    for (double az : azs) {
      for (double el : els) {
        Vec3d uv = azel_to_unitvec(az, el);
        REQUIRE_THAT(uv.norm(), WithinAbs(1.0, 1e-12));
        Vec2d back = unitvec_to_azel(uv);
        REQUIRE_THAT(back[0], WithinAbs(az, 1e-9));  // azimuth
        REQUIRE_THAT(back[1], WithinAbs(el, 1e-9));  // elevation
      }
    }
  }
}

TEST_CASE("environment.plasma.raytrace.compute_ne") {
  CwdGuard guard;
  SetupPlasma();

  RayTraceConfig config = FortranConfig();
  double t = SampleTj2000();

  // Positions are ECEF in kilometers (compute_ne divides by RE[km]). Sample a
  // low F-region point and two higher topside/plasmasphere points on the same ray.
  auto ne_at_alt_km = [&](double alt_km) {
    Vec3d pos(RE + alt_km, 0.0, 0.0);  // on the equatorial x-axis
    return compute_ne(t, pos, config);
  };

  SECTION("electron density is finite, non-negative, and decreases with altitude") {
    double ne_low = ne_at_alt_km(400.0);     // F-region
    double ne_mid = ne_at_alt_km(3000.0);    // topside
    double ne_high = ne_at_alt_km(12000.0);  // plasmasphere

    for (double ne : {ne_low, ne_mid, ne_high}) {
      REQUIRE(std::isfinite(ne));
      REQUIRE(ne >= 0.0);
    }
    // Daytime F-region density is positive.
    REQUIRE(ne_low > 0.0);
    // Total electron density falls off monotonically with geocentric distance.
    REQUIRE(ne_low >= ne_mid);
    REQUIRE(ne_mid >= ne_high);
  }

  SECTION("an unset Kp index is rejected") {
    RayTraceConfig bad = config;
    bad.kp = -1.0;  // default sentinel: Kp not provided
    Vec3d pos(RE + 400.0, 0.0, 0.0);
    REQUIRE_THROWS(compute_ne(t, pos, bad));
  }
}
