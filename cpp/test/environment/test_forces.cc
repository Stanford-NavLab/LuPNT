#include <lupnt/environment/forces.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("environment.forces") {
  SECTION("AccelerationRelativisticCorrection matches circular-orbit form") {
    Real radius = 7000.0e3;
    Real speed = sqrt(GM_EARTH / radius);
    Vec3 r(radius, 0.0, 0.0);
    Vec3 v(0.0, speed, 0.0);

    Vec3 a_newton = -GM_EARTH * r / pow(radius, 3);
    Vec3 expected = a_newton * (3.0 * speed * speed / (C * C));
    Vec3 actual = AccelerationRelativisticCorrection(r, v, GM_EARTH);

    for (int i = 0; i < 3; ++i) {
      REQUIRE_THAT(actual(i).val(), WithinAbs(expected(i).val(), 1.0e-18));
    }
  }

  SECTION("ShadowFunction returns full illumination on the dayside") {
    Vec3 r(R_EARTH + 700.0e3, 0.0, 0.0);
    Vec3 r_sun(AU, 0.0, 0.0);

    REQUIRE_THAT(ShadowFunction(r, r_sun, R_EARTH).val(), WithinAbs(1.0, epsilon));
    REQUIRE_THAT(Illumination(r, r_sun, R_EARTH).val(), WithinAbs(1.0, epsilon));
  }

  SECTION("ShadowFunction returns zero in umbra") {
    Vec3 r(-(R_EARTH + 700.0e3), 0.0, 0.0);
    Vec3 r_sun(AU, 0.0, 0.0);

    REQUIRE_THAT(ShadowFunction(r, r_sun, R_EARTH).val(), WithinAbs(0.0, epsilon));
    REQUIRE_THAT(Illumination(r, r_sun, R_EARTH).val(), WithinAbs(0.0, epsilon));
  }

  SECTION("ShadowFunction returns fractional illumination in penumbra") {
    Vec3 r(-42164.0e3, 6.30e6, 0.0);
    Vec3 r_sun(AU, 0.0, 0.0);

    Real nu = ShadowFunction(r, r_sun, R_EARTH);
    REQUIRE(nu.val() > 0.0);
    REQUIRE(nu.val() < 1.0);
    REQUIRE_THAT(Illumination(r, r_sun, R_EARTH).val(), WithinAbs(nu.val(), epsilon));
  }

  SECTION("ShadowFunction handles maximum partial occultation after the umbra vertex") {
    Vec3 r(-2.0e9, 0.0, 0.0);
    Vec3 r_sun(AU, 0.0, 0.0);

    Real a = asin(R_SUN / (r_sun - r).norm());
    Real b = asin(R_EARTH / r.norm());
    Real expected = 1.0 - b * b / (a * a);

    REQUIRE(b.val() < a.val());
    REQUIRE_THAT(ShadowFunction(r, r_sun, R_EARTH).val(), WithinAbs(expected.val(), epsilon));
  }
}
