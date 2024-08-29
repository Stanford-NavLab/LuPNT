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

  template <typename T> NBodyDynamics<T>::NBodyDynamics(IntegratorType integ)
      : NumericalOrbitDynamics(std::bind(&NBodyDynamics::ComputeRates, this, std::placeholders::_1,
                                         std::placeholders::_2),
                               integ){};
  template class NBodyDynamics<double>;
  template class NBodyDynamics<Real>;

  template <typename T> Vec6 NBodyDynamics<T>::ComputeRates(Real t_tai, const Vec6& rv) const {
    if (frame_ == Frame::NONE) throw std::runtime_error("Frame not set");

    // Position, velocity, and acceleration [km, km/s, km/s^2]
    // w.r.t. to the inertial frame origin
    Vec3 r = rv.head(3);
    Vec3 v = rv.tail(3);
    Vec3 a = Vec3::Zero();

    for (const auto& body : bodies_) {
      if (body.use_gravity_field) {
        auto& grav = body.gravity_field;
        // Position (body-fixed) [km]
        Vector<T, 3> r_bf = ConvertFrame(t_tai, r, frame_, body.fixed_frame).template cast<T>();
        // Acceleration (body-fixed) [km/s^2]
        Vec3 a_bf = AccelarationGravityField<T>(r_bf, grav.GM, grav.R, grav.CS, grav.n, grav.m);
        // Acceleration (inertial) [km/s^2]
        Vec3 ai = ConvertFrame(t_tai, a_bf, body.fixed_frame, frame_, true);
        a += ai;
      } else {
        // Body position w.r.t. the inertial frame origin [km]
        Vec3 r_body = GetBodyPosVel(t_tai, body.id, frame_).head(3);
        // Acceleration (inertial) [km/s^2]
        Vec3 ai = AccelerationPointMass(rv.head(3), r_body, body.GM);
        a += ai;
      }

      // Solar radiation pressure
      if (use_srp_ && body.id != NaifId::SUN) {
        Vec3 r_sun = GetBodyPosVel(t_tai, body.id, NaifId::SUN, frame_).head(3);
        Vec3 a_srp = Illumination(r, r_sun, body.R)
                     * AccelerationSolarRadiation(r, r_sun, area_, mass_, CR_, P_SUN, AU);
        a += a_srp;
      }

      // Atmospheric drag
      if (use_drag_ && body.id == NaifId::EARTH) {
        // TODO: Currently only works for Earth
        Real tt = ConvertTime(t_tai, Time::TAI, Time::TT);
        Real mjd_tt = (tt + MJD_J2000) / SECS_DAY;
        MatX3 Rot = NutationMatrix(mjd_tt) * PrecessionMatrix(MJD_J2000, mjd_tt);
        Vec3 a_drag = AccelerationDrag(mjd_tt, rv, Rot, area_, mass_, CD_);
        a += a_drag;
      }
    }

    Vec6 rv_dot;
    rv_dot << v, a;
    return rv_dot;
  }

  template <typename T> OrbitState NBodyDynamics<T>::PropagateState(const OrbitState& state,
                                                                    Real t0, Real tf, Mat6d* stm) {
    assert(state.GetOrbitStateRepres() == OrbitStateRepres::CARTESIAN
           && "OrbitState type not supported");
    Vec6 xf = Propagate(state.GetVec(), t0, tf, stm);
    return CartesianOrbitState(xf, state.GetFrame());
  }

}  // namespace lupnt
