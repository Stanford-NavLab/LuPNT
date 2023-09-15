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
#include <autodiff/forward/real.hpp>

#include "Dynamics.h"
#include "lupnt/core/Constants.h"
#include "lupnt/dynamics/GravityField.h"
#include "lupnt/physics/SpiceInterface.h"

namespace ad = autodiff;

namespace LPT {

Body Body::Moon(int n_max, int m_max) {
  Body moon;
  moon.name = "MOON";
  moon.id = BodyId::MOON;
  moon.sphericalHarmonics = n_max > 0 || m_max > 0;
  moon.normalized = true;
  moon.n_max = n_max;
  moon.m_max = m_max;

  BodyData bd = GetBodyData(moon.id);
  moon.mu = bd.GM;
  moon.R = bd.R;

  if (moon.sphericalHarmonics)
    std::tie(moon.Cnm, moon.Snm) = LoadGravityCoefficients(bd, n_max);
  return moon;
}

Body Body::Earth(int n_max, int m_max) {
  Body earth;
  earth.name = "EARTH";
  earth.id = BodyId::EARTH;
  earth.sphericalHarmonics = n_max > 0 || m_max > 0;
  earth.normalized = true;
  earth.n_max = n_max;
  earth.m_max = m_max;

  BodyData bd = GetBodyData(earth.id);
  earth.mu = bd.GM;
  earth.R = bd.R;

  if (earth.sphericalHarmonics)
    std::tie(earth.Cnm, earth.Snm) = LoadGravityCoefficients(bd, n_max);
  return earth;
}

NBodyDynamics::NBodyDynamics(std::string integratorType)
    : NumericalDynamics(std::bind(&NBodyDynamics::ComputeRates, this,
                                  std::placeholders::_1, std::placeholders::_2),
                        OrbitStateRepres::CARTESIAN, integratorType){};

void NBodyDynamics::SetCentralBody(const Body &body) { centralBody = body; }

void NBodyDynamics::AddBody(const Body &body) { bodies.push_back(body); }

ad::VectorXreal NBodyDynamics::ComputeRates(ad::real t,
                                            const ad::VectorXreal &x) const {
  ad::VectorXreal acc = ad::VectorXreal::Zero(6);

  // N-body gravity
  acc.head(3) = x.tail(3);
  acc.tail(3) += ComputeNBodyGravity(t, x);

  // Solar radiation pressure
  ad::Vector3real r_body2sc = x.head(3);
  ad::Vector3real r_body2sun =
      SpiceInterface::GetBodyPos("SUN", t, "J2000", "MOON", "NONE");
  ad::Vector3real r_sun2sc = r_body2sc - r_body2sun;
  double R_body = centralBody.R;
  double R_SUN = 696000.0;
  double m = 1000.0;
  double CR = 1.5;
  double area = 1.0;
  ad::Vector3real a_srp = ComputeSolarRadiationPressure(
      r_body2sc, r_sun2sc, R_body, R_SUN, m, CR, area);
  acc.tail(3) += a_srp;

  return acc;
}

ad::Vector3real NBodyDynamics::ComputeNBodyGravity(
    ad::real t, const ad::VectorXreal &rv) const {
  assert(rv.size() == 6);

  ad::Vector3real r = rv.head(3);  // s/c position w.r.t. the center body [km]
  ad::Vector3real v = rv.tail(3);  // s/c velocity w.r.t. the center body [km/s]
  ad::Vector3real a = ad::Vector3real::Zero();  // s/c acceleration w.r.t. the
                                                // center body [km/s^2]

  for (Body body : bodies) {
    ad::VectorXreal rv_i, rv_body;
    ad::Vector3real a_i, a_i_C, r_i;

    // Check central body
    if (body.id != centralBody.id) {
      // i-th body pos and vel w.r.t. the center body [km, km/s]
      rv_body =
          SpiceInterface::GetBodyPosVel(t, (int)centralBody.id, (int)body.id);

      // s/c pos and vel w.r.t. the i-th body [km, km/s]
      rv_i = rv - rv_body;
      r_i = rv_i.head(3);
      a_i_C = -body.mu * rv_body.head(3) / pow((rv_body.head(3)).norm(), 3);  //
    } else {
      rv_i = rv;
      r_i = rv_i.head(3);
      a_i_C = ad::Vector3real::Zero();
    }

    // Check spherical harmonics
    if (body.sphericalHarmonics) {
      std::string inertialFrame = "J2000";
      std::string rotatingFrame = "IAU_" + body.name;
      Eigen::Matrix3d Ur2j = SpiceInterface::GetFrameConversionMatrix(
                                 t, rotatingFrame, inertialFrame)
                                 .block(0, 0, 3, 3);
      ad::Vector3real r_i_rot = Ur2j.transpose() * r_i;
      ad::Vector3real a_i_rot = spharm_acc_ecr(
          body.n_max, body.m_max, r_i_rot, body.R, body.mu, body.Cnm, body.Snm);
      a_i = Ur2j * a_i_rot;
    } else {
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
ad::Vector3real NBodyDynamics::ComputeSolarRadiationPressure(
    const ad::Vector3real &r_body2sc, const ad::Vector3real &r_sun2sc,
    double R_body, double R_SUN, double m, double CR, double area) const {
  // Acceleration of a satellite due to the solar radiation pressure
  ad::Vector3real a_srp;

  // Apparent radius of the occulted body (i.e. the Sun)
  ad::real a = asin(R_SUN / r_sun2sc.norm());
  // Apparent radius of the occulting body
  ad::real b = asin(R_body / r_body2sc.norm());
  // Apparent separation of the centers of both bodies
  ad::real c =
      acos(r_body2sc.dot(r_sun2sc) / (r_body2sc.norm() * r_sun2sc.norm()));

  ad::real nu;
  if (abs(a - b) < c && c < (a + b)) {
    // Occultated area
    ad::real x = (c * c + a * a - b * b) / (2 * c);
    ad::real y = sqrt(a * a - x * x);
    ad::real A = a * a * acos(x / a) + b * b * acos((c - x) / b) - c * y;

    // Remaining fraction of Sun light
    nu = 1 - A / (PI * a * a);
  } else {
    nu = 1;
  }
  a_srp = nu * P_SUN * CR * (area / m) * (r_sun2sc / pow(r_sun2sc.norm(), 3)) *
          (AU * AU);
  return a_srp;
}

}  // namespace LPT