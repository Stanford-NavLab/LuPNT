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

  SECTION("AccelerationPointMass matches the direct-plus-indirect closed form") {
    Vec3 r(7000.0e3, 1500.0e3, -800.0e3);        // spacecraft w.r.t. central body
    Vec3 s(384400.0e3, 120000.0e3, 60000.0e3);   // third body (e.g. Moon) w.r.t. central body
    Real GM = GM_MOON;

    Vec3 d = r - s;
    Vec3 expected = -GM * (d / pow(d.norm(), 3) + s / pow(s.norm(), 3));
    Vec3 actual = AccelerationPointMass(r, s, GM);

    for (int i = 0; i < 3; ++i)
      REQUIRE_THAT(actual(i).val(), WithinAbs(expected(i).val(), 1e-18));
  }

  SECTION("AccelerationSolarRadiation points away from the Sun with inverse-square magnitude") {
    Vec3 r(R_EARTH + 700.0e3, 0.0, 0.0);
    Vec3 r_sun(AU, 0.0, 0.0);
    Real bcoeff = 0.02, P0 = 4.56e-6;

    Vec3 a = AccelerationSolarRadiation(r, r_sun, bcoeff, P0, AU);
    Vec3 d = r - r_sun;  // Sun -> spacecraft; SRP pushes along this direction
    Real expected_mag = bcoeff * P0 * (AU * AU) / d.squaredNorm();

    REQUIRE_THAT(a.norm().val(), WithinRel(expected_mag.val(), 1e-12));
    // Direction: unit(a) == unit(d) (away from the Sun, i.e. -x here).
    Vec3 ua = a / a.norm(), ud = d / d.norm();
    for (int i = 0; i < 3; ++i) REQUIRE_THAT(ua(i).val(), WithinAbs(ud(i).val(), 1e-12));
    REQUIRE(a(0).val() < 0.0);
  }

  SECTION("AccelarationGravityField reduces to point mass at degree/order zero") {
    Vec3 r(7000.0e3, 1000.0e3, -500.0e3);
    Real GM = GM_EARTH, R_ref = R_EARTH;
    MatX CS = MatX::Zero(4, 4);
    CS(0, 0) = 1.0;  // C00 = 1, all higher harmonics zero

    Vec3 actual = AccelarationGravityField<Real>(r, GM, R_ref, CS, 0, 0);
    Vec3 expected = -GM * r / pow(r.norm(), 3);

    for (int i = 0; i < 3; ++i)
      REQUIRE_THAT(actual(i).val(), WithinRel(expected(i).val(), 1e-9));
  }

  SECTION("DensityHarrisPriester is positive and decreases with altitude") {
    Real mjd_tt = 51544.5;  // J2000 epoch
    Vec3 r_low(R_EARTH + 300.0e3, 0.0, 0.0);
    Vec3 r_high(R_EARTH + 600.0e3, 0.0, 0.0);

    Real rho_low = DensityHarrisPriester(mjd_tt, r_low);
    Real rho_high = DensityHarrisPriester(mjd_tt, r_high);

    REQUIRE(rho_low.val() > 0.0);
    REQUIRE(rho_high.val() > 0.0);
    REQUIRE(rho_low.val() > rho_high.val());
  }
}
