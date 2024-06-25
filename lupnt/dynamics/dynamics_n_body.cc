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
#include "lupnt/dynamics/gravity_field.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/body.h"
#include "lupnt/physics/spice_interface.h"

namespace lupnt {

NBodyDynamics::NBodyDynamics(std::string integratorType)
    : NumericalOrbitDynamics(
          std::bind(&NBodyDynamics::ComputeRates, this, std::placeholders::_1,
                    std::placeholders::_2),
          OrbitStateRepres::CARTESIAN, integratorType) {};

VectorX NBodyDynamics::ComputeRates(real t_tai, const VectorX &rv) const {
  VectorX rv_dot = VectorX::Zero(6);

  // N-body gravity
  rv_dot.head(3) = rv.segment(3, 3);  // extract velocity
  rv_dot.tail(3) += ComputeNBodyGravity(t_tai, rv.head(3));

  // Solar radiation pressure
  if (use_srp_) {
    Vector3 r_body2sc = rv.head(3);
    Vector3 r_body2sun = SpiceInterface::GetBodyPosVel(
                             t_tai, NaifId::MOON, NaifId::SUN, Frame::MOON_CI)
                             .head(3);
    Vector3 r_sun2sc = r_body2sc - r_body2sun;
    double R_body = central_body_.R;
    real B_srp = CR_ * (area_ / mass_);  // ballistic coefficient [m^2/kg]
    Vector3 a_srp =
        ComputeSolarRadiationPressure(r_body2sc, r_sun2sc, B_srp, R_body);
    rv_dot.tail(3) += a_srp;
  }

  return rv_dot;
}

Vector3 NBodyDynamics::ComputeNBodyGravity(real t_tai, const Vector3 &r) const {
  assert(r.size() == 3);
  Vector3 a = Vector3::Zero();  // s/c acceleration w.r.t. the
                                // center body [km/s^2]

  for (Body body : bodies_) {
    VectorX rv_i, r_body;
    Vector3 a_i, a_i_C, r_i;

    // Correction for non-central body
    if (body.id != central_body_.id) {
      // i-th body pos and vel w.r.t. the center body [km, km/s]
      r_body = SpiceInterface::GetBodyPosVel(t_tai, central_body_.id, body.id,
                                             Frame::MOON_CI)
                   .head(3);

      // s/c pos and vel w.r.t. the i-th body [km, km/s]
      r_i = r - r_body.head(3);
      a_i_C = -body.mu * r_body / pow((r_body).norm(), 3);  //
    } else {
      r_i = r;
      a_i_C = Vector3::Zero();
    }

    // Check spherical harmonics
    if (body.sphericalHarmonics) {
      // Matrix3d Ur2j = SpiceInterface::GetFrameConversionMatrix(
      //                     t_tai, body.fixed_frame, Frame::GCRF)
      //                     .block(0, 0, 3, 3);
      Matrix3d Ur2j = Matrix3d::Identity();
      Vector3 r_i_rot = Ur2j.transpose() * r_i;
      Vector3 a_i_rot = spharm_acc_ecr(body.n_max, body.m_max, r_i_rot, body.R,
                                       body.mu, body.Cnm, body.Snm);
      a_i = Ur2j * a_i_rot;
    } else {
      // no spherical harmonics
      a_i = -body.mu * r_i / pow((r_i).norm(), 3);
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
Vector3 NBodyDynamics::ComputeSolarRadiationPressure(const Vector3 &r_body2sc,
                                                     const Vector3 &r_sun2sc,
                                                     const real B_srp,
                                                     double R_body) const {
  // Acceleration of a satellite due to the solar radiation pressure
  Vector3 a_srp;

  double R_SUN = 696000.0;  // [km]

  // Apparent radius of the occulted body (i.e. the Sun)
  real a = asin(R_SUN / r_sun2sc.norm());
  // Apparent radius of the occulting body
  real b = asin(R_body / r_body2sc.norm());
  // Apparent separation of the centers of both bodies_
  real c =
      safe_acos(r_body2sc.dot(r_sun2sc) / (r_body2sc.norm() * r_sun2sc.norm()));

  real nu;

  if (abs(a - b) < c && c < (a + b)) {
    // Occultated area
    real x = (c * c + a * a - b * b) / (2 * c);
    real y = sqrt(a * a - x * x);
    real A = a * a * safe_acos(x / a) + b * b * safe_acos((c - x) / b) - c * y;

    // Remaining fraction of Sun light
    nu = 1 - A / (PI * a * a);
  } else {
    nu = 1;
  }
  a_srp = nu * P_SUN * B_srp * (r_sun2sc / pow(r_sun2sc.norm(), 3)) * (AU * AU);
  return a_srp;
}

}  // namespace lupnt