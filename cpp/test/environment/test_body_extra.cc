#include <lupnt/environment/body.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// These tests exercise the spherical-harmonic gravity-field loading path that the
// point-mass-only cases in test_body.cc leave uncovered: the n/m > 1 factory
// overloads (Body::Moon/Earth), which internally call ReadHarmonicGravityField to
// populate GravityField::{GM,R,CS}. (ReadHarmonicGravityField's own template
// instantiation is not exported from the shared library, so it is exercised here
// through the factory constructors that consume it.) Anchors are closed-form: the
// unnormalized zonal C20 = -sqrt(5) * Cbar_20, and J2 = -C20, so both bodies have
// C20 < 0 and a J2 that matches the well-known reference values (Earth ~1.0826e-3,
// Moon ~2.03e-4).

namespace {
  constexpr double kEps = 1e-6;

  // Reference-radius / GM values stored in the bundled coefficient-file headers
  // (these differ slightly from the top-level GM_MOON/R_MOON constants, which is
  // exactly why the field carries its own GM and R).
  constexpr double kMoonFileGM = 4.90279996708864e12;  // [m^3/s^2]
  constexpr double kMoonFileR = 1.738e6;               // [m]
  constexpr double kEarthFileGM = 3.986004415e14;      // [m^3/s^2]
  constexpr double kEarthFileR = 6.3781363e6;          // [m]

  // Unnormalized J2 = -C20 = sqrt(5) * (-Cbar_20).
  double J2FromField(const GravityField<Real>& gf) { return -gf.CS(2, 0).val(); }
}  // namespace

TEST_CASE("environment.body.gravity_field") {
  SECTION("Moon factory loads a GRGM900C spherical-harmonic field") {
    const int n = 8, m = 8;
    Body moon = Body::Moon(n, m, "grgm900c.cof");

    REQUIRE(moon.id == BodyId::MOON);
    REQUIRE(moon.use_gravity_field);

    const GravityField<Real>& gf = moon.gravity_field;
    // Requested truncation is stored, and never exceeds the field's own maximum.
    REQUIRE(gf.n == n);
    REQUIRE(gf.m == m);
    REQUIRE(gf.n_max >= n);
    REQUIRE(gf.m_max >= m);
    // CS is (n+1) x (m+1) with the conventional C00 = 1.
    REQUIRE(gf.CS.rows() == n + 1);
    REQUIRE(gf.CS.cols() == m + 1);
    REQUIRE_THAT(gf.CS(0, 0).val(), WithinAbs(1.0, 1e-12));

    // GM and reference radius come from the file header.
    REQUIRE_THAT(gf.GM.val(), WithinRel(kMoonFileGM, 1e-9));
    REQUIRE_THAT(gf.R.val(), WithinRel(kMoonFileR, 1e-9));
    // ...and are consistent with the canonical constant to <0.1%.
    REQUIRE_THAT(gf.GM.val(), WithinRel(GM_MOON, 1e-3));

    // Oblateness: C20 < 0, and J2 matches the reference Moon value ~2.03e-4.
    REQUIRE(gf.CS(2, 0).val() < 0.0);
    REQUIRE_THAT(J2FromField(gf), WithinRel(2.032e-4, 5e-3));
  }

  SECTION("Earth factory loads an EGM96 spherical-harmonic field") {
    const int n = 8, m = 8;
    Body earth = Body::Earth(n, m, "EGM96.cof");

    REQUIRE(earth.id == BodyId::EARTH);
    REQUIRE(earth.use_gravity_field);

    const GravityField<Real>& gf = earth.gravity_field;
    REQUIRE(gf.CS.rows() == n + 1);
    REQUIRE(gf.CS.cols() == m + 1);
    REQUIRE_THAT(gf.GM.val(), WithinRel(kEarthFileGM, 1e-9));
    REQUIRE_THAT(gf.R.val(), WithinRel(kEarthFileR, 1e-9));
    REQUIRE_THAT(gf.GM.val(), WithinRel(GM_EARTH, 1e-3));

    // J2 for Earth is the textbook 1.0826e-3, and C20 is negative.
    REQUIRE(gf.CS(2, 0).val() < 0.0);
    REQUIRE_THAT(J2FromField(gf), WithinRel(1.08263e-3, 1e-3));
  }

  SECTION("truncation degree/order controls the coefficient-matrix size") {
    GravityField<Real> low = Body::Earth(2, 2, "EGM96.cof").gravity_field;
    GravityField<Real> high = Body::Earth(10, 10, "EGM96.cof").gravity_field;

    REQUIRE(low.CS.rows() == 3);
    REQUIRE(low.CS.cols() == 3);
    REQUIRE(high.CS.rows() == 11);
    REQUIRE(high.CS.cols() == 11);
    // GM/R are field metadata, independent of the requested truncation.
    REQUIRE_THAT(low.GM.val(), WithinRel(high.GM.val(), 1e-12));
    // The shared low-degree coefficients are identical regardless of truncation.
    REQUIRE_THAT(low.CS(2, 0).val(), WithinAbs(high.CS(2, 0).val(), 1e-14));
    REQUIRE_THAT(low.CS(2, 2).val(), WithinAbs(high.CS(2, 2).val(), 1e-14));
  }

  SECTION("degree/order above the field maximum is rejected") {
    // grgm900c is a 360x360 field; asking for more must throw (the exception
    // propagates out of ReadHarmonicGravityField through the factory).
    REQUIRE_THROWS(Body::Moon(400, 400, "grgm900c.cof"));
  }

  SECTION("low degree/order requests stay point-mass") {
    // n or m <= 1 disables the field (use_gravity_field == n_max>1 && m_max>1).
    Body moon0 = Body::Moon(0, 0, "grgm900c.cof");
    Body moon1 = Body::Moon(1, 1, "grgm900c.cof");
    REQUIRE_FALSE(moon0.use_gravity_field);
    REQUIRE_FALSE(moon1.use_gravity_field);
    REQUIRE_THAT(moon0.GM.val(), WithinRel(GM_MOON, kEps));
  }

  SECTION("kilometer units scale the loaded gravity field consistently") {
    Body moon_si = Body::Moon(SI_UNITS, 4, 4, "grgm900c.cof");
    Body moon_km = Body::Moon(KM_S_KG_UNITS, 4, 4, "grgm900c.cof");

    REQUIRE(moon_km.use_gravity_field);
    REQUIRE(moon_km.gravity_field.units == KM_S_KG_UNITS);
    // GM scales as length^3 (m^3 -> km^3 = /1e9); R scales as length (/1e3).
    REQUIRE_THAT(moon_km.gravity_field.GM.val(),
                 WithinRel(moon_si.gravity_field.GM.val() / 1.0e9, kEps));
    REQUIRE_THAT(moon_km.gravity_field.R.val(),
                 WithinRel(moon_si.gravity_field.R.val() / 1.0e3, kEps));
    // Dimensionless coefficients are unit-independent.
    REQUIRE_THAT(moon_km.gravity_field.CS(2, 0).val(),
                 WithinAbs(moon_si.gravity_field.CS(2, 0).val(), 1e-14));
  }
}
