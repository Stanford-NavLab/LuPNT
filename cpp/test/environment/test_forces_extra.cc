#include <lupnt/environment/forces.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  const double eps_extra = 1e-6;
}

// ---------------------------------------------------------------------------
// AccelarationGravityField: the J2 (degree-2 zonal) contribution must match the
// closed-form oblateness perturbation of Montenbruck & Gill. Uses the standard
// unnormalized-coefficient convention C_{2,0} = -J2 (and C_{0,0}=1 for the
// central point mass), so the total should equal (point mass + closed-form J2).
// ---------------------------------------------------------------------------
TEST_CASE("environment.forces_extra.j2_zonal") {
  const Real J2 = 1.08262668e-3;  // Earth's second zonal harmonic
  const Real GM = GM_EARTH, R = R_EARTH;

  MatX CS = MatX::Zero(3, 3);
  CS(0, 0) = 1.0;  // central point mass
  CS(2, 0) = -J2;  // C_{2,0} = -J2 (unnormalized)

  auto j2_closed_form = [&](const Vec3& r) {
    Real rn = r.norm();
    Vec3 a_pm = -GM * r / pow(rn, 3);
    Real f = -1.5 * J2 * (GM / (rn * rn)) * pow(R / rn, 2);
    Real z2r2 = pow(r(2) / rn, 2);
    Vec3 a_j2(f * (1.0 - 5.0 * z2r2) * r(0) / rn, f * (1.0 - 5.0 * z2r2) * r(1) / rn,
              f * (3.0 - 5.0 * z2r2) * r(2) / rn);
    return (a_pm + a_j2).eval();
  };

  SECTION("general position matches the closed-form J2 acceleration") {
    Vec3 r(7000.0e3, 1200.0e3, 3500.0e3);
    Vec3 actual = AccelarationGravityField<Real>(r, GM, R, CS, 2, 0);
    Vec3 expected = j2_closed_form(r);
    for (int i = 0; i < 3; ++i) REQUIRE_THAT(actual(i).val(), WithinRel(expected(i).val(), 1e-6));
  }

  SECTION("equatorial point: J2 perturbation is purely radial and inward") {
    Vec3 r(7000.0e3, 0.0, 0.0);  // on the +x equator
    Vec3 a_full = AccelarationGravityField<Real>(r, GM, R, CS, 2, 0);
    Vec3 a_pm = -GM * r / pow(r.norm(), 3);
    Vec3 dJ2 = a_full - a_pm;  // the J2 contribution
    // The extra pull from the equatorial bulge is inward (toward -x) on the
    // equator, and lies along the radial (x) axis only.
    REQUIRE(dJ2(0).val() < 0.0);
    REQUIRE_THAT(dJ2(1).val(), WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(dJ2(2).val(), WithinAbs(0.0, 1e-12));
  }

  SECTION("north-pole point: J2 perturbation pushes outward along +z") {
    Vec3 r(0.0, 0.0, 7000.0e3);  // over the north pole
    Vec3 a_full = AccelarationGravityField<Real>(r, GM, R, CS, 2, 0);
    Vec3 a_pm = -GM * r / pow(r.norm(), 3);
    Vec3 dJ2 = a_full - a_pm;
    REQUIRE(dJ2(2).val() > 0.0);
    REQUIRE_THAT(dJ2(0).val(), WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(dJ2(1).val(), WithinAbs(0.0, 1e-12));
  }
}

// ---------------------------------------------------------------------------
// AccelerationDrag: with an identity inertial->body transform the acceleration
// must oppose the atmosphere-relative velocity (v - omega x r) and scale
// linearly with the ballistic coefficient. Physical invariants only -- no
// dependence on the absolute (unit-laden) drag magnitude.
// ---------------------------------------------------------------------------
TEST_CASE("environment.forces_extra.drag") {
  Real mjd_tt = 51544.5;  // J2000 epoch
  // Prograde near-circular LEO state (equatorial), so v_rel is well-defined.
  Real radius = R_EARTH + 400.0e3;
  Real speed = sqrt(GM_EARTH / radius);
  Vec6 rv;
  rv << radius, 0.0, 0.0, 0.0, speed, 0.0;
  Mat3 T = Mat3::Identity();

  const Vec3 omega(0.0, 0.0, 7.29212e-5);
  Vec3 v_rel = rv.tail(3) - omega.cross(Vec3(rv.head(3)));

  SECTION("drag acceleration is anti-parallel to the atmosphere-relative velocity") {
    Real bcoeff = 0.02;
    Vec3 a = AccelerationDrag(mjd_tt, rv, T, bcoeff);

    // Opposes relative velocity: a . v_rel < 0.
    REQUIRE(a.dot(v_rel).val() < 0.0);
    // Exactly anti-parallel: a x v_rel == 0.
    Vec3 cross = Vec3(a).cross(v_rel);
    for (int i = 0; i < 3; ++i) REQUIRE_THAT(cross(i).val(), WithinAbs(0.0, 1e-20));
  }

  SECTION("drag acceleration scales linearly with the ballistic coefficient") {
    Vec3 a1 = AccelerationDrag(mjd_tt, rv, T, Real(0.02));
    Vec3 a2 = AccelerationDrag(mjd_tt, rv, T, Real(0.04));
    for (int i = 0; i < 3; ++i) REQUIRE_THAT(a2(i).val(), WithinRel(2.0 * a1(i).val(), 1e-12));
  }

  SECTION("zero ballistic coefficient yields zero drag") {
    Vec3 a = AccelerationDrag(mjd_tt, rv, T, Real(0.0));
    for (int i = 0; i < 3; ++i) REQUIRE_THAT(a(i).val(), WithinAbs(0.0, 1e-30));
  }
}

// ---------------------------------------------------------------------------
// AccelerationEarthSpacecraft: the composite Earth-spacecraft acceleration
// (central harmonic gravity + luni-solar point masses + relativity + SRP +
// drag) is dominated by the central gravity term, whose magnitude is exactly
// mu/r^2 for a degree-0 field (the body-frame rotation preserves its norm).
// The remaining perturbations are ~1e-5 m/s^2 or smaller, so the total
// magnitude must stay within a small relative band of mu/r^2. Uses only
// low-precision analytic Sun/Moon positions (no ephemeris/SPICE required).
// ---------------------------------------------------------------------------
TEST_CASE("environment.forces_extra.earth_spacecraft_gravity_dominated") {
  Real mjd_tt = 51544.5;  // J2000
  Real radius = R_EARTH + 800.0e3;
  Real speed = sqrt(GM_EARTH / radius);
  Vec6 rv;
  rv << radius, 0.0, 0.0, 0.0, speed, 0.0;

  // Degree/order 0 gravity field == central point mass.
  GravityField<Real> grav;
  grav.n_max = 0;
  grav.m_max = 0;
  grav.n = 0;
  grav.m = 0;
  grav.GM = GM_EARTH;
  grav.R = R_EARTH;
  grav.CS = MatX::Zero(1, 1);
  grav.CS(0, 0) = 1.0;

  SECTION("total magnitude is dominated by central gravity mu/r^2") {
    Vec3 a = AccelerationEarthSpacecraft(mjd_tt, rv, Real(0.0), Real(0.0), grav);
    Real expected_mag = GM_EARTH / (radius * radius);
    // Perturbations (luni-solar ~1e-5, relativity ~1e-8 m/s^2) are >5 orders of
    // magnitude below the ~8.6 m/s^2 central term, so |a| stays within ~1e-5.
    REQUIRE_THAT(a.norm().val(), WithinRel(expected_mag.val(), 1e-4));
    // Finite result.
    for (int i = 0; i < 3; ++i) REQUIRE(std::isfinite(a(i).val()));
  }

  SECTION("enabling SRP/drag perturbs the total but keeps it gravity-dominated") {
    Vec3 a_bare = AccelerationEarthSpacecraft(mjd_tt, rv, Real(0.0), Real(0.0), grav);
    Vec3 a_full = AccelerationEarthSpacecraft(mjd_tt, rv, Real(0.02), Real(0.02), grav);
    // SRP+drag change the acceleration...
    REQUIRE((a_full - a_bare).norm().val() > 0.0);
    // ...but only slightly relative to central gravity.
    Real expected_mag = GM_EARTH / (radius * radius);
    REQUIRE((a_full - a_bare).norm().val() < 1e-3 * expected_mag.val());
  }
}
