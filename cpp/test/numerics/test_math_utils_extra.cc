#include <lupnt/numerics/math_utils.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {

  // RatioOfSectorToTriangleArea(r1, r2, tau) returns Gauss's sector-to-triangle
  // area ratio eta, with tau = sqrt(GM) * dt (Montenbruck & Eberhard convention,
  // see CartToClassical in state_conversions.cc).
  //
  // Closed-form anchor: on a *circular* orbit of radius R, the area swept along
  // the orbit between two radii separated by central angle theta is the circular
  // sector 1/2 R^2 theta, while the straight-line triangle they form with the
  // focus has area 1/2 R^2 sin(theta). Hence eta = theta / sin(theta) exactly,
  // independent of R and GM.
  TEST_CASE("numerics.math_utils_extra.sector_triangle_ratio_circular") {
    const double R = 7000.0e3;
    const double GM = GM_EARTH;
    const double n = std::sqrt(GM / (R * R * R));  // circular mean motion

    auto eta_for_theta = [&](double theta) {
      Vec3 r1(R, 0.0, 0.0);
      Vec3 r2(R * std::cos(theta), R * std::sin(theta), 0.0);
      double dt = theta / n;          // time of flight along circular orbit
      Real tau = std::sqrt(GM) * dt;  // Montenbruck's scaled time argument
      return RatioOfSectorToTriangleArea(r1, r2, tau).val();
    };

    // theta / sin(theta) closed form, checked at several arc lengths.
    for (double theta : {0.2, 0.5, 1.0, 1.5}) {
      double expected = theta / std::sin(theta);
      REQUIRE_THAT(eta_for_theta(theta), WithinAbs(expected, 1e-6));
    }

    // eta is strictly > 1 for any non-degenerate arc, and increases with arc
    // length (the sector bulges further beyond the chord triangle).
    REQUIRE(eta_for_theta(1.0) > eta_for_theta(0.5));
    REQUIRE(eta_for_theta(0.5) > 1.0);
  }

  TEST_CASE("numerics.math_utils_extra.sector_triangle_ratio_small_arc") {
    // As the arc shrinks to zero, sector and triangle areas coincide: eta -> 1.
    const double R = 1737.4e3;
    const double GM = GM_MOON;
    const double n = std::sqrt(GM / (R * R * R));

    double theta = 1e-3;
    Vec3 r1(R, 0.0, 0.0);
    Vec3 r2(R * std::cos(theta), R * std::sin(theta), 0.0);
    Real tau = std::sqrt(GM) * (theta / n);
    Real eta = RatioOfSectorToTriangleArea(r1, r2, tau);
    REQUIRE_THAT(eta.val(), WithinAbs(1.0, 1e-6));
  }

  // Independence from the common length scale: eta depends only on the geometry
  // (theta) and the time-of-flight consistency, not on R. Re-derive at a very
  // different radius and confirm the same theta/sin(theta) value.
  TEST_CASE("numerics.math_utils_extra.sector_triangle_ratio_scale_invariance") {
    const double GM = GM_EARTH;
    const double theta = 0.8;
    const double expected = theta / std::sin(theta);

    for (double R : {6.6e6, 4.2e7}) {  // LEO-ish and GEO-ish radii
      const double n = std::sqrt(GM / (R * R * R));
      Vec3 r1(R, 0.0, 0.0);
      Vec3 r2(R * std::cos(theta), R * std::sin(theta), 0.0);
      Real tau = std::sqrt(GM) * (theta / n);
      REQUIRE_THAT(RatioOfSectorToTriangleArea(r1, r2, tau).val(), WithinAbs(expected, 1e-6));
    }
  }

}  // namespace
