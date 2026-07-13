#include <lupnt/environment/solar_system.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// The existing physics/test_solar_system.cc *calls* every solar-system helper but
// asserts nothing on the results. This file adds physical-invariant checks:
// rotation matrices must be orthogonal with unit determinant, inertial<->body-fixed
// rotations must round-trip to the identity, the pos-vel (6x6) rotations must carry
// the 3x3 attitude in their diagonal blocks, and the analytic Sun/Moon positions
// must sit at physically sane geocentric distances.

namespace {
  // Verify M is a proper rotation: M M^T = I and det(M) = +1.
  void RequireOrthonormal(const Mat3& M, double tol = 1e-9) {
    Mat3 I = M * M.transpose();
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) REQUIRE_THAT(I(i, j).val(), WithinAbs(i == j ? 1.0 : 0.0, tol));
    REQUIRE_THAT(M.determinant().val(), WithinAbs(1.0, tol));
  }

  // A representative TDB epoch (~2025) expressed as seconds since J2000, plus the
  // matching MJD(TT) used by the obliquity/nutation helpers.
  const Real kTdb = 25.0 * 365.25 * SECS_DAY;  // ~2025 in seconds past J2000
  const Real kMjdTt = MJD_J2000_TT + 25.0 * 365.25;
}  // namespace

TEST_CASE("environment.solar_system.rotation_invariants") {
  SECTION("nutation matrices are proper rotations and nearly agree") {
    Mat3 R_nut = NutationMatrix(kMjdTt);
    Mat3 R_nut_low = NutationMatrixLowPrecision(kMjdTt);
    RequireOrthonormal(R_nut);
    RequireOrthonormal(R_nut_low);
    // Full and low-precision nutation differ only at the sub-arcsecond level,
    // so every element matches to well under 1e-4.
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        REQUIRE_THAT(R_nut(i, j).val(), WithinAbs(R_nut_low(i, j).val(), 1e-4));
  }

  SECTION("nutation angles are small and finite") {
    auto [dpsi, deps] = NutAngles(kMjdTt);
    REQUIRE(std::isfinite(dpsi.val()));
    REQUIRE(std::isfinite(deps.val()));
    // Nutation in longitude/obliquity is at most ~20 arcsec.
    REQUIRE(std::abs(dpsi.val()) < 30.0 * RAD_ARCSEC);
    REQUIRE(std::abs(deps.val()) < 30.0 * RAD_ARCSEC);
    // Equation of the equinoxes = dpsi * cos(eps) is consistent in magnitude.
    Real eqeq = EquinoxEquation(kMjdTt);
    REQUIRE(std::abs(eqeq.val()) <= std::abs(dpsi.val()) + 1e-12);
  }

  SECTION("equatorial-to-ecliptic and Greenwich rotations are proper") {
    Mat3 R_eq2ecl = Equatorial2EclipticMatrix(kMjdTt);
    Mat3 R_gha = GreenwichHourAngleMatrix(kMjdTt);
    RequireOrthonormal(R_eq2ecl);
    RequireOrthonormal(R_gha);
    // Equatorial->ecliptic is a pure rotation about the X axis by the obliquity,
    // so the X row/column is preserved.
    Real eps = MeanObliquity(kMjdTt);
    REQUIRE_THAT(R_eq2ecl(0, 0).val(), WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(R_eq2ecl(1, 1).val(), WithinAbs(std::cos(eps.val()), 1e-12));
    // Greenwich hour-angle rotation is about Z, preserving the Z axis.
    REQUIRE_THAT(R_gha(2, 2).val(), WithinAbs(1.0, 1e-12));
  }

  SECTION("precession over a nonzero interval is a proper rotation") {
    Mat3 P = PrecessionMatrix(MJD_J2000_TT, kMjdTt);
    RequireOrthonormal(P);
    // 25 years of precession is a small angle: the matrix stays close to identity.
    REQUIRE_THAT(P(0, 0).val(), WithinAbs(1.0, 1e-3));
  }
}

TEST_CASE("environment.solar_system.body_frame_rotations") {
  const std::vector<BodyId> bodies = {BodyId::MARS, BodyId::JUPITER, BodyId::VENUS};

  SECTION("inertial<->body-fixed rotations are inverses") {
    for (BodyId id : bodies) {
      Mat3 R_i2b = RotPosInertialToBodyFixed(id, kTdb);
      Mat3 R_b2i = RotPosBodyFixedToInertial(id, kTdb);
      RequireOrthonormal(R_i2b);
      RequireOrthonormal(R_b2i);
      // Body-fixed<-inertial is the transpose/inverse of inertial<-body-fixed.
      Mat3 prod = R_i2b * R_b2i;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          REQUIRE_THAT(prod(i, j).val(), WithinAbs(i == j ? 1.0 : 0.0, 1e-9));
    }
  }

  SECTION("pos-vel rotations carry the attitude in their diagonal block") {
    for (BodyId id : bodies) {
      Mat3 R_i2b = RotPosInertialToBodyFixed(id, kTdb);
      Mat6 R6_i2b = RotPosVelInertialToBodyFixed(id, kTdb);
      Mat6 R6_b2i = RotPosVelBodyFixedToInertial(id, kTdb);

      // Upper-left 3x3 (position block) equals the pure-position rotation.
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          REQUIRE_THAT(R6_i2b(i, j).val(), WithinAbs(R_i2b(i, j).val(), 1e-9));

      // The lower-right 3x3 (velocity attitude) matches the position block.
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          REQUIRE_THAT(R6_i2b(3 + i, 3 + j).val(), WithinAbs(R_i2b(i, j).val(), 1e-9));

      // The 6x6 inertial<->body-fixed transforms are mutual inverses.
      Mat6 prod = R6_i2b * R6_b2i;
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
          REQUIRE_THAT(prod(i, j).val(), WithinAbs(i == j ? 1.0 : 0.0, 1e-6));
    }
  }

  SECTION("planet orientation angles are finite") {
    for (BodyId id : bodies) {
      Vec4 a = PlanetOrientation(id, kTdb);
      for (int i = 0; i < 4; ++i) REQUIRE(std::isfinite(a(i).val()));
    }
  }
}

TEST_CASE("environment.solar_system.analytic_body_geometry") {
  SECTION("low-precision Moon distance is within the lunar orbit range") {
    // Geocentric Moon distance stays inside [356500, 406700] km physically;
    // use a slightly wider band to cover the low-precision series.
    Vec3 r_moon = MoonPositionLowPrecision(kMjdTt);
    double d = r_moon.norm().val();
    REQUIRE(d > 3.5e8);
    REQUIRE(d < 4.1e8);
  }

  SECTION("low-precision Sun distance is about one AU") {
    Vec3 r_sun = SunPositionLowPrecision(kMjdTt);
    double d = r_sun.norm().val();
    // Earth-Sun distance varies ~+-1.7% about 1 AU over the year.
    REQUIRE_THAT(d, WithinRel(AU, 0.03));
  }
}
