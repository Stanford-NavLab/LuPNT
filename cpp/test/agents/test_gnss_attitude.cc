#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 2e-8;

TEST_CASE("agents.gnss_attitude") {
  // Representative GPS satellite & Sun positions, in ECI [m]
  Vec3 r_sat(26560e3, 0.0, 0.0);
  Vec3 r_sun(0.3 * AU, 0.9 * AU, 0.1 * AU);

  Vec3 ex_ref, ey_ref, ez_ref;
  GnssAttitude::Compute(r_sat, r_sun, ex_ref, ey_ref, ez_ref);

  SECTION("Compute produces an orthonormal, right-handed, nadir-pointing triad") {
    // Each axis is a unit vector
    RequireNear(Real(ex_ref.norm()), Real(1.0), epsilon);
    RequireNear(Real(ey_ref.norm()), Real(1.0), epsilon);
    RequireNear(Real(ez_ref.norm()), Real(1.0), epsilon);

    // Mutually orthogonal
    RequireNear(ex_ref.dot(ey_ref), Real(0.0), epsilon);
    RequireNear(ey_ref.dot(ez_ref), Real(0.0), epsilon);
    RequireNear(ez_ref.dot(ex_ref), Real(0.0), epsilon);

    // Right-handed: ex = ey x ez
    RequireNear(ex_ref, ey_ref.cross(ez_ref), epsilon);

    // ez points to nadir, i.e. along -r_sat
    RequireNear(ez_ref, (-r_sat).normalized(), epsilon);
  }

  SECTION("Constructor immediately computes the same triad as the static Compute") {
    GnssAttitude att(r_sat, r_sun);
    RequireNear(att.GetEx(), ex_ref, epsilon);
    RequireNear(att.GetEy(), ey_ref, epsilon);
    RequireNear(att.GetEz(), ez_ref, epsilon);
  }

  SECTION("Default-constructed instance starts from the canonical axes") {
    GnssAttitude att;
    RequireNear(att.GetEx(), Vec3(Vec3::UnitX()), epsilon);
    RequireNear(att.GetEy(), Vec3(Vec3::UnitY()), epsilon);
    RequireNear(att.GetEz(), Vec3(Vec3::UnitZ()), epsilon);
  }

  SECTION("Update (re-)computes and overwrites the cached triad") {
    GnssAttitude att;
    att.Update(r_sat, r_sun);
    RequireNear(att.GetEx(), ex_ref, epsilon);
    RequireNear(att.GetEy(), ey_ref, epsilon);
    RequireNear(att.GetEz(), ez_ref, epsilon);

    // A different geometry produces a different (and still valid) triad
    Vec3 r_sat2(0.0, 26560e3, 5000e3);
    Vec3 r_sun2(-0.5 * AU, 0.2 * AU, 0.7 * AU);
    Vec3 ex2, ey2, ez2;
    GnssAttitude::Compute(r_sat2, r_sun2, ex2, ey2, ez2);

    att.Update(r_sat2, r_sun2);
    RequireNear(att.GetEx(), ex2, epsilon);
    RequireNear(att.GetEy(), ey2, epsilon);
    RequireNear(att.GetEz(), ez2, epsilon);
  }

  SECTION("GetRotationMatrix assembles [ex, ey, ez] into an orthonormal rotation matrix") {
    GnssAttitude att(r_sat, r_sun);
    Mat3 R = att.GetRotationMatrix();

    RequireNear(Vec3(R.col(0)), ex_ref, epsilon);
    RequireNear(Vec3(R.col(1)), ey_ref, epsilon);
    RequireNear(Vec3(R.col(2)), ez_ref, epsilon);

    // Orthonormal: R^T * R = I
    Mat3 RtR = R.transpose() * R;
    Mat3 I = Mat3::Identity();
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) RequireNear(RtR(i, j), I(i, j), epsilon);
    }

    // Right-handed (proper rotation): det(R) = +1
    RequireNear(Real(R.determinant()), Real(1.0), epsilon);
  }

  SECTION("GetAngles returns the boresight-relative (theta, phi) of a unit direction") {
    GnssAttitude att(r_sat, r_sun);
    Real theta, phi;

    // Looking straight down boresight (along ez): off-boresight angle phi = 0
    att.GetAngles(att.GetEz(), theta, phi);
    RequireNear(phi, Real(0.0), epsilon);

    // Looking along ex: on the phi = 90 deg cone, at azimuth theta = 0
    att.GetAngles(att.GetEx(), theta, phi);
    RequireNear(phi, Real(PI_OVER_TWO), epsilon);
    RequireNear(theta, Real(0.0), epsilon);

    // Looking along ey: phi = 90 deg, azimuth theta = 90 deg
    att.GetAngles(att.GetEy(), theta, phi);
    RequireNear(phi, Real(PI_OVER_TWO), epsilon);
    RequireNear(theta, Real(PI_OVER_TWO), epsilon);

    // Looking opposite boresight (along -ez): phi = 180 deg
    att.GetAngles(Vec3(-att.GetEz()), theta, phi);
    RequireNear(phi, Real(PI), epsilon);
  }
}

TEST_CASE("agents.gnss_attitude_yaw_steering") {
  // Representative GPS satellite state (position + velocity, ECI [m], [m/s])
  // and Sun position (ECI [m]). `v_sat` is roughly circular-orbit speed for a
  // GPS-altitude orbit, chosen perpendicular to `r_sat` so (r_sat, v_sat) span
  // a sensible orbital plane.
  Vec3 r_sat(26560e3, 0.0, 0.0);
  Vec3 v_sat(0.0, 3873.7, 0.0);
  Vec3 r_sun(0.3 * AU, 0.9 * AU, 0.1 * AU);

  // Looser tolerance for cross-checking two different floating-point paths
  // (geometric Sun-pointing vs. trig-chain through GnssYawSteering) that are
  // mathematically identical but numerically distinct.
  const double cross_path_epsilon = 1e-7;

  SECTION(
      "Yaw-law-driven Compute(r_sat, v_sat, r_sun, ...) matches the "
      "Sun-direction-based Compute(r_sat, r_sun, ...)") {
    // This is the key connection verified numerically while planning this
    // change: the existing geometric Sun-pointing frame is exactly the
    // orbital reference frame rotated about nadir by the nominal yaw angle
    // phi_nom = GnssYawSteering::NominalYawAngle(beta, mu).
    Vec3 ex_geom, ey_geom, ez_geom;
    GnssAttitude::Compute(r_sat, r_sun, ex_geom, ey_geom, ez_geom);

    Vec3 ex_yaw, ey_yaw, ez_yaw;
    GnssAttitude::Compute(r_sat, v_sat, r_sun, ex_yaw, ey_yaw, ez_yaw);

    RequireNear(ex_yaw, ex_geom, cross_path_epsilon);
    RequireNear(ey_yaw, ey_geom, cross_path_epsilon);
    RequireNear(ez_yaw, ez_geom, cross_path_epsilon);
  }

  SECTION("ComputeFromYawAngle with the nominal yaw angle reproduces the Sun-pointing frame") {
    Vec3 ex_geom, ey_geom, ez_geom;
    GnssAttitude::Compute(r_sat, r_sun, ex_geom, ey_geom, ez_geom);

    Real beta = GnssYawSteering::BetaAngle(r_sat, v_sat, r_sun);
    Real mu = GnssYawSteering::OrbitAngle(r_sat, v_sat, r_sun);
    Real phi_nom = GnssYawSteering::NominalYawAngle(beta, mu);

    Vec3 ex_phi, ey_phi, ez_phi;
    GnssAttitude::ComputeFromYawAngle(r_sat, v_sat, phi_nom, ex_phi, ey_phi, ez_phi);

    RequireNear(ex_phi, ex_geom, cross_path_epsilon);
    RequireNear(ey_phi, ey_geom, cross_path_epsilon);
    RequireNear(ez_phi, ez_geom, cross_path_epsilon);
  }

  SECTION(
      "ComputeFromYawAngle produces an orthonormal, right-handed, nadir-pointing triad "
      "for an arbitrary yaw angle") {
    Real phi = Real(0.37);  // arbitrary yaw angle [rad], not the nominal value
    Vec3 ex, ey, ez;
    GnssAttitude::ComputeFromYawAngle(r_sat, v_sat, phi, ex, ey, ez);

    RequireNear(Real(ex.norm()), Real(1.0), epsilon);
    RequireNear(Real(ey.norm()), Real(1.0), epsilon);
    RequireNear(Real(ez.norm()), Real(1.0), epsilon);

    RequireNear(ex.dot(ey), Real(0.0), epsilon);
    RequireNear(ey.dot(ez), Real(0.0), epsilon);
    RequireNear(ez.dot(ex), Real(0.0), epsilon);

    RequireNear(ex, ey.cross(ez), epsilon);
    RequireNear(ez, (-r_sat).normalized(), epsilon);
  }

  SECTION("Update(r_sat, v_sat, r_sun) caches the same triad as the yaw-law-driven Compute") {
    Vec3 ex_ref, ey_ref, ez_ref;
    GnssAttitude::Compute(r_sat, v_sat, r_sun, ex_ref, ey_ref, ez_ref);

    GnssAttitude att;
    att.Update(r_sat, v_sat, r_sun);
    RequireNear(att.GetEx(), ex_ref, epsilon);
    RequireNear(att.GetEy(), ey_ref, epsilon);
    RequireNear(att.GetEz(), ez_ref, epsilon);
  }
}
