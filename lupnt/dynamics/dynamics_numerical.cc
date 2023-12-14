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

#include "dynamics.h"
#include "lupnt/core/constants.h"

namespace lupnt {

// ****************************************************************************
// NumericalDynamics
// ****************************************************************************

void NumericalDynamics::Propagate(OrbitState &state, real t0, real tf,
                                  real dt) {
  assert(state.GetOrbitStateRepres() == state_representation_);
  Vector6 x0 = state.GetVector();
  Vector6 xf = propagator_.Propagate(odefunc_, t0, tf, x0, dt);
  state.SetVector(xf);
}

void NumericalDynamics::Propagate(Vector6 &x, real t0, real tf, real dt) {
  Vector6 x0 = x;
  Vector6 xf = propagator_.Propagate(odefunc_, t0, tf, x0, dt);
  x = xf;
}

void NumericalDynamics::PropagateWithStm(OrbitState &state, real t0, real tf,
                                         real dt, Matrix6d &stm) {
  assert(state.GetOrbitStateRepres() == state_representation_);
  Vector6 x0 = state.GetVector();
  MatrixXd J(6, 6);
  VectorX xf = propagator_.PropagateWithStm(odefunc_, t0, tf, x0, dt, J);
  state.SetVector(xf);
  stm = J;
}

void NumericalDynamics::PropagateWithStm(Vector6 &x, real t0, real tf, real dt,
                                         Matrix6d &stm) {
  Vector6 x0 = x;
  MatrixXd J(6, 6);
  VectorX xf = propagator_.PropagateWithStm(odefunc_, t0, tf, x0, dt, J);
  x = xf;
  stm = J;
}

void NumericalDynamics::Propagate(OrbitState &state, real t0, real tf) {
  real dt_prop = tf - t0;
  if (dt_ == 0.0) {
    dt_prop = (tf - t0) / 10;
  } else {
    dt_prop = dt_;
  }
  Propagate(state, t0, tf, dt_prop);
}

void NumericalDynamics::Propagate(Vector6 &x, real t0, real tf) {
  real dt_prop = tf - t0;
  if (dt_ == 0.0) {
    dt_prop = (tf - t0) / 10;
  } else {
    dt_prop = dt_;
  }
  Propagate(x, t0, tf, dt_prop);
}

void NumericalDynamics::PropagateWithStm(OrbitState &state, real t0, real tf,
                                         Matrix6d &stm) {
  real dt_prop = tf - t0;
  if (dt_ == 0.0) {
    dt_prop = (tf - t0) / 10;
  } else {
    dt_prop = dt_;
  }
  PropagateWithStm(state, t0, tf, dt_prop, stm);
}

void NumericalDynamics::PropagateWithStm(Vector6 &x, real t0, real tf,
                                         Matrix6d &stm) {
  real dt_prop = tf - t0;
  if (dt_ == 0.0) {
    dt_prop = (tf - t0) / 10;
  } else {
    dt_prop = dt_;
  }
  PropagateWithStm(x, t0, tf, dt_prop, stm);
}

// arbitrary state size
void NumericalDynamics::PropagateX(VectorX &x, real t0, real tf) {
  real dt_prop = tf - t0;
  if (dt_ == 0.0) {
    dt_prop = (tf - t0) / 10;
  } else {
    dt_prop = dt_;
  }

  Vector6 x6;
  x6 << x(0), x(1), x(2), x(3), x(4), x(5);
  Propagate(x6, t0, tf, dt_prop);
  x.head(6) = x6;
}

void NumericalDynamics::PropagateWithStmX(VectorX &x, real t0, real tf,
                                          MatrixXd &stm) {
  real dt_prop = tf - t0;
  if (dt_ == 0.0) {
    dt_prop = (tf - t0) / 10;
  } else {
    dt_prop = dt_;
  }

  Vector6 x6;
  x6 << x(0), x(1), x(2), x(3), x(4), x(5);
  Matrix6d stm6;
  stm6 = stm.block(0, 0, 6, 6);
  PropagateWithStm(x6, t0, tf, dt_prop, stm6);
  x.head(6) = x6;
  stm.block(0, 0, 6, 6) = stm6;
}

// ****************************************************************************
// CartesianTwoBodyDynamics
// ****************************************************************************

CartesianTwoBodyDynamics::CartesianTwoBodyDynamics(double mu,
                                                   std::string integratorType)
    : mu_(mu),
      NumericalDynamics(std::bind(&CartesianTwoBodyDynamics::ComputeRates, this,
                                  std::placeholders::_1, std::placeholders::_2),
                        OrbitStateRepres::CARTESIAN, integratorType){};

VectorX CartesianTwoBodyDynamics::ComputeRates(real t, const VectorX &x) const {
  Vector6 dxdt;
  Vector3 r = x.head(3);
  Vector3 v = x.tail(3);
  real r_norm = r.norm();

  dxdt.head(3) = v;
  dxdt.tail(3) = -mu_ * r / pow(r_norm, 3);

  return dxdt;
}

// ****************************************************************************
// MoonFixedDynamics
// ****************************************************************************

MoonFixedDynamics::MoonFixedDynamics(double mu, std::string integratorType)
    : mu_(mu),
      NumericalDynamics(std::bind(&MoonFixedDynamics::ComputeRates, this,
                                  std::placeholders::_1, std::placeholders::_2),
                        OrbitStateRepres::CARTESIAN, integratorType){};

VectorX MoonFixedDynamics::ComputeRates(real t, const VectorX &x) const {
  Vector6 dxdt;
  Vector3 r = x.head(3);
  Vector3 v = x.tail(3);
  real r_norm = r.norm();

  dxdt.head(3) = Vector3::Zero();
  dxdt.tail(3) = Vector3::Zero();

  return dxdt;
}

// ****************************************************************************
// J2CartesianTwoBodyDynamics
// ****************************************************************************

J2CartesianTwoBodyDynamics::J2CartesianTwoBodyDynamics(
    double mu, double J2, double Rbody, std::string integratorType)
    : mu_(mu),
      J2_(J2),
      Rbody_(Rbody),
      NumericalDynamics(
          std::bind(&J2CartesianTwoBodyDynamics::ComputeRates, this,
                    std::placeholders::_1, std::placeholders::_2),
          OrbitStateRepres::CARTESIAN, integratorType){};

VectorX J2CartesianTwoBodyDynamics::ComputeRates(real t,
                                                 const VectorX &x) const {
  VectorX acc(6);

  Vector3 r = x.head(3);
  Vector3 v = x.tail(3);
  real r_norm = r.norm();

  acc.head(3) = v;
  acc.tail(3) = -mu_ * r / pow(r_norm, 3);

  Vector3 aJ2;
  real aux1 = -3.0 / 2.0 * mu_ * J2_ * pow(Rbody_, 2.0) / pow(r_norm, 4.0);
  real aux2 = 5.0 * pow(r(2) / Rbody_, 2.0);
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

VectorX J2KeplerianDynamics::ComputeRates(real t, const VectorX &x) const {
  real p = x(0) * (1.0 - x(1) * x(1));
  real n = sqrt(mu / pow(x(0), 3.0));
  real eta = sqrt(1.0 - x(1) * x(1));

  VectorX coeDot(6);
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

VectorX MoonMeanDynamics::ComputeRates(real t, const VectorX &x) const {
  real a = x(0);
  real e = x(1);
  real i = x(2);
  real O = x(3);
  real w = x(4);
  real M = x(5);

  VectorX coeDot(6);
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
};  // namespace lupnt
