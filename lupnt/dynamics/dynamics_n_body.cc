/**
 * @file NBodyDynamics.cpp
 * @author Stanford NAV LAB
 * @brief Multiple-body dynamics
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "dynamics.h"
#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/body.h"
#include "lupnt/physics/spice_interface.h"

namespace lupnt {

NBodyDynamics::NBodyDynamics(std::string integratorType)
    : NumericalOrbitDynamics(
          std::bind(&NBodyDynamics::ComputeRates, this, std::placeholders::_1,
                    std::placeholders::_2),
          OrbitStateRepres::CARTESIAN, integratorType) {};

VecX NBodyDynamics::ComputeRates(Real t_tai, const VecX &rv) const {
  VecX rv_dot = VecX::Zero(6);

  // N-body gravity
  rv_dot.head(3) = rv.segment(3, 3);  // extract velocity
  rv_dot.tail(3) += ComputeNBodyGravity(t_tai, rv.head(3));

  // Solar radiation pressure
  if (use_srp_) {
    Vec3 r_body2sc = rv.head(3);
    Vec3 r_body2sun =
        GetBodyPosVel(t_tai, NaifId::MOON, NaifId::SUN, Frame::MOON_CI).head(3);
    Vec3 r_sun2sc = r_body2sc - r_body2sun;
    double R_body = central_body_.R;
    Real B_srp = CR_ * (area_ / mass_);  // ballistic coefficient [m^2/kg]
    Vec3 a_srp =
        ComputeSolarRadiationPressure(r_body2sc, r_sun2sc, B_srp, R_body);
    rv_dot.tail(3) += a_srp;
  }

  return rv_dot;
}

Vec3 NBodyDynamics::ComputeNBodyGravity(Real t_tai, const Vec3 &r) const {
  assert(r.size() == 3);
  Vec3 a = Vec3::Zero();  // s/c acceleration w.r.t. the
                          // center body [km/s^2]

  for (Body body : bodies_) {
    VecX rv_i, r_body;
    Vec3 a_i, a_i_C, r_i;

    // Correction for non-central body
    if (body.id != central_body_.id) {
      // i-th body pos and vel w.r.t. the center body [km, km/s]
      r_body = GetBodyPosVel(t_tai, central_body_.id, body.id, Frame::MOON_CI)
                   .head(3);

      // s/c pos and vel w.r.t. the i-th body [km, km/s]
      r_i = r - r_body.head(3);
      a_i_C = -body.GM * r_body / pow((r_body).norm(), 3);  //
    } else {
      r_i = r;
      a_i_C = Vec3::Zero();
    }

    // Check if body.gravity_field has been set
    if (body.has_gravity_field) {
      // Mat3d Ur2j = GetFrameConversionMat(
      //                     t_tai, body.fixed_frame, Frame::GCRF)
      //                     .block(0, 0, 3, 3);
      Mat3d Ur2j = Mat3d::Identity();
      Vec3 r_i_rot = Ur2j.transpose() * r_i;
      Vec3 a_i_rot = Vec3::Zero();
      a_i = Ur2j * a_i_rot;
    } else {
      // no spherical harmonics
      a_i = -body.GM * r_i / pow((r_i).norm(), 3);
    }

    a += a_i + a_i_C;
  }

  return a;
}

/**
 * @brief
 *
 * @ref 1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and
 * Applications_, 2012, p.77-83.
 */
Vec3 NBodyDynamics::ComputeSolarRadiationPressure(const Vec3 &r_body2sc,
                                                  const Vec3 &r_sun2sc,
                                                  const Real B_srp,
                                                  double R_body) const {
  // Acceleration of a satellite due to the solar radiation pressure
  Vec3 a_srp;

  double R_SUN = 696000.0;  // [km]

  // Apparent radius of the occulted body (i.e. the Sun)
  Real a = asin(R_SUN / r_sun2sc.norm());
  // Apparent radius of the occulting body
  Real b = asin(R_body / r_body2sc.norm());
  // Apparent separation of the centers of both bodies_
  Real c =
      safe_acos(r_body2sc.dot(r_sun2sc) / (r_body2sc.norm() * r_sun2sc.norm()));

  Real nu;

  if (abs(a - b) < c && c < (a + b)) {
    // Occultated area
    Real x = (c * c + a * a - b * b) / (2 * c);
    Real y = sqrt(a * a - x * x);
    Real A = a * a * safe_acos(x / a) + b * b * safe_acos((c - x) / b) - c * y;

    // Remaining fraction of Sun light
    nu = 1 - A / (PI * a * a);
  } else {
    nu = 1;
  }
  a_srp = nu * P_SUN * B_srp * (r_sun2sc / pow(r_sun2sc.norm(), 3)) * (AU * AU);
  return a_srp;
}

}  // namespace lupnt