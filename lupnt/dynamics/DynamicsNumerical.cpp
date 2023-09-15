/**
 * @file DynamicsNumerical.cpp
 * @author Stanford NAV LAB
 * @brief List of Numerical Dynamics
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "Dynamics.h"
#include "lupnt/core/Constants.h"

namespace LPT {

// ****************************************************************************
// NumericalDynamics
// ****************************************************************************

void NumericalDynamics::Propagate(OrbitState &state, ad::real t0, ad::real tf,
                                  ad::real dt) {
  assert(state.GetOrbitStateRepres() == stateRepres_);
  ad::Vector6real x0 = state.GetVector();
  ad::Vector6real xf = propagator_.Propagate(odefunc_, t0, tf, x0, dt);
  state.SetVector(xf);
}

void NumericalDynamics::Propagate(ad::Vector6real &x, ad::real t0, ad::real tf,
                                  ad::real dt) {
  ad::Vector6real x0 = x;
  ad::Vector6real xf = propagator_.Propagate(odefunc_, t0, tf, x0, dt);
  x = xf;
}

void NumericalDynamics::PropagateWithStm(OrbitState &state, ad::real t0,
                                         ad::real tf, ad::real dt,
                                         Eigen::Matrix6d &stm) {
  assert(state.GetOrbitStateRepres() == stateRepres_);
  ad::Vector6real x0 = state.GetVector();
  Eigen::MatrixXd J(6, 6);
  ad::VectorXreal xf =
      propagator_.PropagateWithStm(odefunc_, t0, tf, x0, dt, J);
  state.SetVector(xf);
  stm = J;
}

void NumericalDynamics::PropagateWithStm(ad::Vector6real &x, ad::real t0,
                                         ad::real tf, ad::real dt,
                                         Eigen::Matrix6d &stm) {
  ad::Vector6real x0 = x;
  Eigen::MatrixXd J(6, 6);
  ad::VectorXreal xf =
      propagator_.PropagateWithStm(odefunc_, t0, tf, x0, dt, J);
  x = xf;
  stm = J;
}

void NumericalDynamics::Propagate(OrbitState &state, ad::real t0, ad::real tf) {
  Propagate(state, t0, tf, dt_);
}

void NumericalDynamics::Propagate(ad::Vector6real &x, ad::real t0,
                                  ad::real tf) {
  Propagate(x, t0, tf, dt_);
}

void NumericalDynamics::PropagateWithStm(OrbitState &state, ad::real t0,
                                         ad::real tf, Eigen::Matrix6d &stm) {
  PropagateWithStm(state, t0, tf, dt_, stm);
}

void NumericalDynamics::PropagateWithStm(ad::Vector6real &x, ad::real t0,
                                         ad::real tf, Eigen::Matrix6d &stm) {
  PropagateWithStm(x, t0, tf, dt_, stm);
}

// ****************************************************************************
// CartesianTwoBodyDynamics
// ****************************************************************************

CartesianTwoBodyDynamics::CartesianTwoBodyDynamics(double mu_in,
                                                   std::string integratorType)
    : mu(mu_in),
      NumericalDynamics(std::bind(&CartesianTwoBodyDynamics::ComputeRates, this,
                                  std::placeholders::_1, std::placeholders::_2),
                        OrbitStateRepres::CARTESIAN, integratorType){};

ad::VectorXreal CartesianTwoBodyDynamics::ComputeRates(
    ad::real t, const ad::VectorXreal &x) const {
  ad::Vector6real dxdt;
  ad::Vector3real r = x.head(3);
  ad::Vector3real v = x.tail(3);
  ad::real r_norm = r.norm();

  dxdt.head(3) = v;
  dxdt.tail(3) = -mu * r / pow(r_norm, 3);

  return dxdt;
}

// ****************************************************************************
// MoonFixedDynamics
// ****************************************************************************

MoonFixedDynamics::MoonFixedDynamics(double mu_in, std::string integratorType)
    : mu(mu_in),
      NumericalDynamics(std::bind(&MoonFixedDynamics::ComputeRates, this,
                                  std::placeholders::_1, std::placeholders::_2),
                        OrbitStateRepres::CARTESIAN, integratorType){};

ad::VectorXreal MoonFixedDynamics::ComputeRates(
    ad::real t, const ad::VectorXreal &x) const {
  ad::Vector6real dxdt;
  ad::Vector3real r = x.head(3);
  ad::Vector3real v = x.tail(3);
  ad::real r_norm = r.norm();

  dxdt.head(3) = ad::Vector3real::Zero();
  dxdt.tail(3) = ad::Vector3real::Zero();

  return dxdt;
}

// ****************************************************************************
// J2CartesianTwoBodyDynamics
// ****************************************************************************

J2CartesianTwoBodyDynamics::J2CartesianTwoBodyDynamics(
    double mu_in, double J2_in, double Rbody_in, std::string integratorType)
    : mu(mu_in),
      J2(J2_in),
      Rbody(Rbody_in),
      NumericalDynamics(
          std::bind(&J2CartesianTwoBodyDynamics::ComputeRates, this,
                    std::placeholders::_1, std::placeholders::_2),
          OrbitStateRepres::CARTESIAN, integratorType){};

ad::VectorXreal J2CartesianTwoBodyDynamics::ComputeRates(
    ad::real t, const ad::VectorXreal &x) const {
  ad::VectorXreal acc(6);

  ad::Vector3real r = x.head(3);
  ad::Vector3real v = x.tail(3);
  ad::real r_norm = r.norm();

  acc.head(3) = v;
  acc.tail(3) = -mu * r / pow(r_norm, 3);

  ad::Vector3real aJ2;
  ad::real aux1 = -3.0 / 2.0 * mu * J2 * Rbody * Rbody / pow(r_norm, 4.0);
  ad::real aux2 = 5.0 * pow(r(2) / Rbody, 2.0);
  aJ2(0) = aux1 * (1.0 - aux2) * r(0);
  aJ2(1) = aux1 * (1.0 - aux2) * r(1);
  aJ2(2) = aux1 * (3.0 - aux2) * r(2);

  acc.tail(3) += aJ2;

  return acc;
}

// ****************************************************************************
// J2KeplerianDynamics
// ****************************************************************************

J2KeplerianDynamics::J2KeplerianDynamics(double mu_in, double J2_in,
                                         double Rbody_in,
                                         std::string integratorType)
    : mu(mu_in),
      J2(J2_in),
      Rbody(Rbody_in),
      NumericalDynamics(std::bind(&J2KeplerianDynamics::ComputeRates, this,
                                  std::placeholders::_1, std::placeholders::_2),
                        OrbitStateRepres::CARTESIAN, integratorType){};

ad::VectorXreal J2KeplerianDynamics::ComputeRates(
    ad::real t, const ad::VectorXreal &x) const {
  ad::real p = x(0) * (1.0 - x(1) * x(1));
  ad::real n = sqrt(mu / pow(x(0), 3.0));
  ad::real eta = sqrt(1.0 - x(1) * x(1));

  ad::VectorXreal coeDot(6);
  coeDot(0) = 0.0;
  coeDot(1) = 0.0;
  coeDot(2) = 0.0;
  coeDot(3) = -3.0 / 2.0 * J2 * pow(Rbody / p, 2.0) * n * cos(x(2));
  coeDot(4) = 3.0 / 4.0 * J2 * pow(Rbody / p, 2.0) * n *
              (5.0 * pow(cos(x(2)), 2.0) - 1.0);
  coeDot(5) = n + 3.0 / 4.0 * J2 * pow(Rbody / p, 2.0) * n * eta *
                      (3.0 * pow(cos(x(2)), 2.0) - 1.0);
  return coeDot;
}

// ****************************************************************************
// MoonMeanDynamics
// ****************************************************************************

MoonMeanDynamics::MoonMeanDynamics(std::string integratorType)
    : NumericalDynamics(std::bind(&MoonMeanDynamics::ComputeRates, this,
                                  std::placeholders::_1, std::placeholders::_2),
                        OrbitStateRepres::CARTESIAN, integratorType){};

ad::VectorXreal MoonMeanDynamics::ComputeRates(ad::real t,
                                               const ad::VectorXreal &x) const {
  ad::real a = x(0);
  ad::real e = x(1);
  ad::real i = x(2);
  ad::real O = x(3);
  ad::real w = x(4);
  ad::real M = x(5);

  ad::VectorXreal coeDot(6);
  coeDot(0) = 0;
  coeDot(1) = (15 * k * pow(n3, 2) * pow(a, 3.0 / 2.0)) / (8 * sqrt(MU_MOON)) *
              e * sqrt(1 - pow(e, 2)) * pow(sin(i), 2) * sin(2 * w);
  coeDot(2) = -(15 * k * pow(n3, 2) * pow(a, 3.0 / 2.0)) /
              (16 * sqrt(MU_MOON)) * pow(e, 2) / sqrt(1 - pow(e, 2)) *
              sin(2 * i) * sin(2 * w);
  coeDot(3) = -(3 * J2 * sqrt(MU_MOON) * pow(R_MOON, 2)) /
                  (2 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 2)) * cos(i) +
              (3 * k * pow(n3, 2) * pow(a, 3.0 / 2.0)) /
                  (8 * sqrt(MU_MOON) * sqrt(1 - pow(e, 2))) *
                  (5 * pow(e, 2) * cos(2 * w) - 3 * pow(e, 2) - 2) * cos(i);
  coeDot(4) = (3 * J2 * sqrt(MU_MOON) * pow(R_MOON, 2)) /
                  (4 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 2)) *
                  (5 * pow(cos(i), 2) - 1) +
              (3 * k * pow(n3, 2) * pow(a, 3.0 / 2.0)) /
                  (8 * sqrt(MU_MOON) * sqrt(1 - pow(e, 2))) *
                  ((5 * pow(cos(i), 2) - 1 + pow(e, 2)) +
                   5 * (1 - pow(e, 2) - pow(cos(i), 2)) * cos(2 * w));
  coeDot(5) = sqrt(MU_MOON) / pow(a, 3.0 / 2.0) +
              (3 * J2 * sqrt(MU_MOON) * pow(R_MOON, 2)) /
                  (4 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 3.0 / 2.0)) *
                  (3 * pow(cos(i), 2) - 1) -
              (k * pow(n3, 2) * pow(a, 3.0 / 2.0)) / (8 * sqrt(MU_MOON)) *
                  ((3 * pow(e, 2) + 7) * (3 * pow(cos(i), 2) - 1) +
                   15 * (1 + pow(e, 2)) * pow(sin(i), 2) * cos(2 * w));
  return coeDot;
}
};  // namespace LPT
