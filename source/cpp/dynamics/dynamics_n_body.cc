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

#include "lupnt/core/constants.h"
#include "lupnt/data/kernels.h"
#include "lupnt/dynamics/dynamics.h"
#include "lupnt/dynamics/forces.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/body.h"
#include "lupnt/physics/solar_system.h"

namespace lupnt {

  NBodyDynamics::NBodyDynamics(std::string integratorType)
      : NumericalOrbitDynamics(std::bind(&NBodyDynamics::ComputeRates, this, std::placeholders::_1,
                                         std::placeholders::_2),
                               OrbitStateRepres::CARTESIAN, integratorType) {};

  VecX NBodyDynamics::ComputeRates(Real t_tai, const VecX& rv) const {
    assert(rv.size() == 6);
    // Position, velocity, and acceleration [km, km/s, km/s^2]
    // w.r.t. to the inertial frame origin
    Vec3 r = rv.head(3);
    Vec3 v = rv.tail(3);
    Vec3 a = Vec3::Zero();

    for (const auto& body : bodies_) {
      if (body.use_gravity_field) {
        auto& grav = body.gravity_field;
        // Position (body-fixed) [km]
        Vec3 r_bf = ConvertFrame(t_tai, r, central_body_.inertial_frame, body.fixed_frame);
        // Acceleration (body-fixed) [km/s^2]
        Vec3 a_bf
            = AccelarationGravityField(r_bf, grav.GM, grav.R, grav.CS, grav.n_max, grav.m_max);
        // Acceleration (inertial) [km/s^2]
        Vec3 ai = ConvertFrame(t_tai, a_bf, body.fixed_frame, central_body_.inertial_frame);
        a += ai;
      } else {
        // Body position w.r.t. the inertial frame origin [km]
        Vec3 r_body = GetBodyPosVel(t_tai, central_body_.id, body.id).head(3);
        // Acceleration (inertial) [km/s^2]
        Vec3 ai = AccelerationPointMass(rv.head(3), r_body, body.GM);
        a += ai;
      }
    }

    // Solar radiation pressure
    if (use_srp_) {
      Vec3 r_sun = GetBodyPosVel(t_tai, central_body_.id, NaifId::SUN).head(3);

      Vec3 a_srp = Illumination(r, r_sun, central_body_.R)
                   * AccelerationSolarRadiation(r, r_sun, area_, mass_, CR_, P_SUN, AU);
      a += a_srp;
    }

    // Atmospheric drag
    if (use_drag_ && central_body_.id == NaifId::EARTH) {
      // TODO: Currently only works for Earth
      Real tt = ConvertTime(t_tai, TimeSys::TAI, TimeSys::TT);
      Real mjd_tt = (tt + MJD_J2000) / SECS_DAY;
      Mat3 T = NutationMatrix(mjd_tt) * PrecessionMatrix(MJD_J2000, mjd_tt);
      Vec3 a_drag = AccelerationDrag(mjd_tt, rv, T, area_, mass_, CD_);
      a += a_drag;
    }

    Vec6 rv_dot;
    rv_dot << v, a;
    return rv_dot;
  }

}  // namespace lupnt
