/**
 * @file gnss_attitude.cc
 * @author Stanford NAV LAB
 * @brief GNSS satellite attitude (body) frame computation
 * @version 0.1
 * @date 2025-06-07
 *
 * @copyright Copyright (c) 2025
 */
#include "lupnt/agents/gnss_attitude.h"

#include "lupnt/agents/gnss_yaw_steering.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {

  GnssAttitude::GnssAttitude(const Vec3& r_sat_eci, const Vec3& r_sun_eci) {
    Update(r_sat_eci, r_sun_eci);
  }

  void GnssAttitude::Compute(const Vec3& r_sat_eci, const Vec3& r_sun_eci, Vec3& ex, Vec3& ey,
                             Vec3& ez) {
    // GNSS satellite attitude frame (nadir / sun cross-track / along-track),
    // mirrors the attitude computation in `GNSSMeas.setup_measurements`:
    //   ez = normalize(-r_sat)                                (nadir)
    //   ey = normalize(cross(ez, normalize(r_sun - r_sat)))   (cross-track, toward Sun)
    //   ex = normalize(cross(ey, ez))                         (along-track, completes triad)
    ez = (-r_sat_eci).normalized();
    Vec3 e_to_sun = (r_sun_eci - r_sat_eci).normalized();
    ey = ez.cross(e_to_sun).normalized();
    ex = ey.cross(ez).normalized();
  }

  void GnssAttitude::ComputeFromYawAngle(const Vec3& r_sat_eci, const Vec3& v_sat_eci,
                                         Real yaw_angle, Vec3& ex, Vec3& ey, Vec3& ez) {
    // Orbital reference frame (the body frame at yaw angle phi = 0):
    //   ez     = normalize(-r_sat)              (nadir)
    //   ex_ref = normalize(v_sat)               (along-track)
    //   ey_ref = cross(ez, ex_ref)              (completes the right-handed triad)
    // The actual frame (ex, ey) is (ex_ref, ey_ref) rotated by `yaw_angle`
    // about the nadir axis ez:
    //   ex =  cos(phi) * ex_ref + sin(phi) * ey_ref
    //   ey = -sin(phi) * ex_ref + cos(phi) * ey_ref
    ez = (-r_sat_eci).normalized();
    Vec3 ex_ref = v_sat_eci.normalized();
    Vec3 ey_ref = ez.cross(ex_ref);
    Real c = cos(yaw_angle);
    Real s = sin(yaw_angle);
    ex = (c * ex_ref + s * ey_ref).normalized();
    ey = (-s * ex_ref + c * ey_ref).normalized();
  }

  void GnssAttitude::Compute(const Vec3& r_sat_eci, const Vec3& v_sat_eci, const Vec3& r_sun_eci,
                             Vec3& ex, Vec3& ey, Vec3& ez) {
    // Drive the attitude frame from the documented nominal yaw-steering law
    // (Kouba, 2009; Eq. 1 of Cheng et al., 2025): compute the Sun-elevation
    // angle `beta` and orbit angle `mu`, evaluate the nominal yaw angle
    // `phi = NominalYawAngle(beta, mu)`, and build the frame by rotating the
    // orbital reference frame about nadir by `phi` (see `ComputeFromYawAngle`).
    // This reproduces exactly the Sun-pointing frame of the geometric
    // `Compute(r_sat_eci, r_sun_eci, ...)` overload above (the nominal yaw law
    // is, by construction, the yaw angle that keeps the solar panels
    // Sun-pointing) -- verified numerically in `test_gnss_attitude.cc`.
    Real beta = GnssYawSteering::BetaAngle(r_sat_eci, v_sat_eci, r_sun_eci);
    Real mu = GnssYawSteering::OrbitAngle(r_sat_eci, v_sat_eci, r_sun_eci);
    Real phi = GnssYawSteering::NominalYawAngle(beta, mu);
    ComputeFromYawAngle(r_sat_eci, v_sat_eci, phi, ex, ey, ez);
  }

  void GnssAttitude::Update(const Vec3& r_sat_eci, const Vec3& r_sun_eci) {
    Compute(r_sat_eci, r_sun_eci, ex_, ey_, ez_);
  }

  void GnssAttitude::Update(const Vec3& r_sat_eci, const Vec3& v_sat_eci, const Vec3& r_sun_eci) {
    Compute(r_sat_eci, v_sat_eci, r_sun_eci, ex_, ey_, ez_);
  }

  Mat3 GnssAttitude::GetRotationMatrix() const {
    Mat3 R;
    R.col(0) = ex_;
    R.col(1) = ey_;
    R.col(2) = ez_;
    return R;
  }

  void GnssAttitude::GetAngles(const Vec3& u, Real& theta, Real& phi) const {
    phi = safe_acos(u.dot(ez_));
    theta = atan2(u.dot(ey_), u.dot(ex_));
  }

}  // namespace lupnt
