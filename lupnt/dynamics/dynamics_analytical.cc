/**
 * @file DynamicsAnalytical.cpp
 * @author Stanford NAV LAB
 * @brief Analytical orbit dynamics implementation
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include <lupnt/physics/orbit_state_utils.h>

#include <Eigen/QR>

#include "dynamics.h"

namespace lupnt {

// ****************************************************************************
// KeplerianDynamics
// ****************************************************************************

KeplerianDynamics::KeplerianDynamics(double mu) : mu_(mu){};

// ClassicalOE
void KeplerianDynamics::Propagate(ClassicalOE &state, real dt) {
  real a = state.a();
  real n = sqrt(mu_ / pow(a, 3));
  state.Set_M(wrapToPi(state.M() + n * dt));
}
void KeplerianDynamics::PropagateWithStm(ClassicalOE &state, real dt,
                                         Matrix6d &stm) {
  real a = state.a();
  real n = sqrt(mu_ / pow(a, 3));
  state.Set_M(wrapToPi(state.M() + n * dt));
  stm = Matrix6d::Identity(6, 6);
  stm(5, 0) = -3.0 / 2.0 * (n / a * dt).val();
}

// QuasiNonsingularOE
void KeplerianDynamics::Propagate(QuasiNonsingularOE &state, real dt) {
  state.Set_u(state.u() + sqrt(mu_ / pow(state.a(), 3)) * dt);
}
void KeplerianDynamics::PropagateWithStm(QuasiNonsingularOE &state, real dt,
                                         Matrix6d &stm) {
  throw std::runtime_error("Not implemented");
}

// EquinoctialOE
void KeplerianDynamics::Propagate(EquinoctialOE &state, real dt) {
  state.Set_lon(state.lon() + sqrt(mu_ / pow(state.a(), 3)) * dt);
}
void KeplerianDynamics::PropagateWithStm(EquinoctialOE &state, real dt,
                                         Matrix6d &stm) {
  throw std::runtime_error("Not implemented");
}

/* ****************************************************************************
  ClohessyWiltshireDynamics
  ************************************************************************** */

ClohessyWiltshireDynamics::ClohessyWiltshireDynamics(real a_in, real n_in)
    : a(a_in), n(n_in){};

void ClohessyWiltshireDynamics::Propagate(OrbitState &state, real tEnd) {
  if (state.GetOrbitStateRepres() != OrbitStateRepres::CARTESIAN)
    throw std::runtime_error("OrbitState type not supported");

  VectorX xEnd = ComputeMatrix(tEnd) * K;
  state.SetVector(xEnd);
}

void ClohessyWiltshireDynamics::Initialize(CartesianOrbitState &state,
                                           real tStart) {
  tInit = tStart;
  MatrixX Phi = ComputeMatrix(tStart);
  K = Phi.colPivHouseholderQr().solve(state.GetVector());
}

MatrixX ClohessyWiltshireDynamics::ComputeMatrix(real t) {
  real sin_nt = sin(n * t);
  real cos_nt = cos(n * t);

  MatrixX A = MatrixX::Zero(6, 6);
  A.block(0, 0, 3, 3) = a * MatrixXd::Identity(3, 3);
  A.block(3, 3, 3, 3) = a * n * MatrixXd::Identity(3, 3);

  MatrixX B = MatrixX::Zero(6, 6);
  B(0, 0) = 1.0;
  B(0, 1) = sin_nt;
  B(0, 2) = cos_nt;

  B(1, 0) = -3.0 / 2.0 * n * t;
  B(1, 1) = 2.0 * cos_nt;
  B(1, 2) = -2.0 * sin_nt;
  B(1, 3) = 1.0;

  B(2, 4) = sin_nt;
  B(2, 5) = cos_nt;

  B(3, 1) = cos_nt;
  B(3, 2) = -sin_nt;

  B(4, 0) = -3.0 / 2.0;
  B(4, 1) = -2.0 * sin_nt;
  B(4, 2) = -2.0 * cos_nt;

  B(5, 4) = cos_nt;
  B(5, 5) = -sin_nt;

  MatrixX Phi = A * B;
  return Phi;
}

// ****************************************************************************
// YamanakaAnkersenDynamics
// ****************************************************************************

YamanakaAnkersenDynamics::YamanakaAnkersenDynamics()
    : a(0.0), n(0.0), e(0.0), M0(0.0){};
void YamanakaAnkersenDynamics::Propagate(CartesianOrbitState &state,
                                         real tEnd) {
  if (state.GetOrbitStateRepres() == OrbitStateRepres::CARTESIAN) {
    VectorX xEnd = ComputeMatrix(tEnd) * K;
    state.SetVector(xEnd);
  } else {
    throw std::runtime_error("OrbitState type not supported");
  }
}
void YamanakaAnkersenDynamics::Initialize(ClassicalOE &coe_c,
                                          CartesianOrbitState &rv_rtn,
                                          real tStart, double mu) {
  a = coe_c.a().val();
  n = sqrt(mu / pow(a, 3.0));
  e = coe_c.e().val();
  M0 = coe_c.M().val();
  tInit = tStart;

  MatrixX Phi = ComputeMatrix(tStart);
  K = ComputeInverseMatrix(tStart) * rv_rtn.GetVector();
  // K = Phi.colPivHouseholderQr().solve(state.GetVector());
}
MatrixX YamanakaAnkersenDynamics::ComputeMatrix(real t) {
  real M = n * (t - tInit) + M0;
  real f = MeanToTrueAnomaly(M, e);
  real sin_f = sin(f);
  real cos_f = cos(f);
  real k = 1.0 + e * cos(f);
  real kp = -e * sin(f);
  real eta = sqrt(1.0 - e * e);
  real tau = n * t / pow(eta, 3.0);

  MatrixX A = MatrixX::Zero(6, 6);
  A.block(0, 0, 3, 3) = a * eta * eta * MatrixXd::Identity(3, 3);
  A.block(3, 3, 3, 3) = a * n / eta * MatrixXd::Identity(3, 3);

  MatrixX B = MatrixX::Zero(6, 6);
  B(0, 0) = 1.0 / k + 3.0 / 2.0 * kp * tau;
  B(0, 1) = sin_f;
  B(0, 2) = cos_f;

  B(1, 0) = -3.0 / 2.0 * k * tau;
  B(1, 1) = (1.0 + 1.0 / k) * cos_f;
  B(1, 2) = -(1.0 + 1.0 / k) * sin_f;
  B(1, 3) = 1.0 / k;

  B(2, 4) = 1.0 / k * sin_f;
  B(2, 5) = 1.0 / k * cos_f;

  B(3, 0) = kp / 2.0 - 3.0 / 2.0 * k * k * (k - 1) * tau;
  B(3, 1) = k * k * cos_f;
  B(3, 2) = -k * k * sin_f;

  B(4, 0) = -3.0 / 2.0 * (k + k * k * kp * tau);
  B(4, 1) = -(k * k + 1) * sin_f;
  B(4, 2) = -e - (k * k + 1) * cos_f;
  B(4, 3) = -kp;

  B(5, 4) = e + cos_f;
  B(5, 5) = -sin_f;

  MatrixX Phi = A * B;
  return Phi;
}
MatrixX YamanakaAnkersenDynamics::ComputeInverseMatrix(real t) {
  real M = n * (t - tInit) + M0;
  real f = MeanToTrueAnomaly(M, e);
  real sin_f = sin(f);
  real cos_f = cos(f);
  real k = 1.0 + e * cos(f);
  real kp = -e * sin(f);
  real eta = sqrt(1.0 - e * e);
  real tau = n * t / pow(eta, 3.0);

  MatrixX A = MatrixX::Zero(6, 6);
  A.block(0, 0, 3, 3) = 1.0 / a / (eta * eta) * MatrixXd::Identity(3, 3);
  A.block(3, 3, 3, 3) = eta / a / n * MatrixXd::Identity(3, 3);

  MatrixX B = MatrixX::Zero(6, 6);
  B(0, 0) = 2.0 * (k * k) * (k + 1) / (eta * eta);
  B(0, 1) = 2.0 * (k * k) * kp / (eta * eta);
  B(0, 3) = -2.0 * kp / (eta * eta);
  B(0, 4) = 2.0 * k / (eta * eta);

  B(1, 0) = (1.0 - pow(k + 1, 2) / (eta * eta)) * sin_f +
            3.0 * e * (k * k) * (k + 1) * tau / (eta * eta);
  B(1, 1) = -(k + 1) * kp * sin_f / (eta * eta) +
            3.0 * e * (k * k) * kp * tau / (eta * eta);
  B(1, 3) = (1.0 / (eta * eta)) * (cos_f - 2.0 * e / k) -
            3.0 * e * kp * tau / (eta * eta);
  B(1, 4) = -(1.0 / (eta * eta)) * (1.0 + 1.0 / k) * sin_f +
            3.0 * e * k * tau / (eta * eta);

  B(2, 0) = -(k / (eta * eta)) * (2.0 * e + (k + 2) * cos_f);
  B(2, 1) = -(1.0 / (eta * eta)) * (e + (k + 1) * cos_f) * kp;
  B(2, 3) = -(1.0 / (eta * eta)) * sin_f;
  B(2, 4) = -(1.0 / (eta * eta)) * (e / k + (1.0 + 1.0 / k) * cos_f);

  B(3, 0) = (pow(k + 1, 2) / (eta * eta)) * kp +
            3.0 * (k * k) * (k + 1) * tau / (eta * eta);
  B(3, 1) = (k / (eta * eta)) * (2.0 + k - (k * k)) - 1.0 +
            3.0 * (k * k) * kp * tau / (eta * eta);
  B(3, 3) =
      (1.0 / (eta * eta)) * (k - 1.0 - 2.0 / k) - 3.0 * kp * tau / (eta * eta);
  B(3, 4) =
      (1.0 / (eta * eta)) * (1.0 + 1.0 / k) * kp + 3.0 * k * tau / (eta * eta);

  B(4, 2) = sin_f;
  B(4, 5) = (1.0 / k) * cos_f;

  B(5, 2) = e + cos_f;
  B(5, 5) = -(1.0 / k) * sin_f;

  MatrixX Phi = B * A;
  return Phi;
}

// ****************************************************************************
// RoeGeometricMappingDynamics
// ****************************************************************************

RoeGeometricMappingDynamics::RoeGeometricMappingDynamics()
    : a(0.0), n(0.0), e(0.0), M0(0.0), ex(0.0), ey(0.0), tInit(0.0){};
void RoeGeometricMappingDynamics::Propagate(CartesianOrbitState &state,
                                            real tEnd) {
  if (state.GetOrbitStateRepres() == OrbitStateRepres::CARTESIAN) {
    std::cout << "xStart: " << state.GetVector().transpose() << std::endl;
    VectorX xEnd = ComputeMatrix(tEnd) * K;
    state.SetVector(xEnd);
    std::cout << "K: " << K.transpose() << std::endl;
    std::cout << "xEnd: " << xEnd.transpose() << std::endl;
  } else {
    throw std::runtime_error("OrbitState type not supported");
  }
}
void RoeGeometricMappingDynamics::Initialize(ClassicalOE coe_c,
                                             QuasiNonsingularROE &roe,
                                             real tStart, double mu) {
  a = coe_c.a().val();
  e = coe_c.e().val();
  i = coe_c.i().val();
  w = coe_c.w().val();
  M0 = coe_c.M().val();

  ex = (double)(coe_c.e() * cos(coe_c.w()));
  ey = (double)(coe_c.e() * sin(coe_c.w()));

  n = (double)sqrt(mu / pow(coe_c.a(), 3));
  K = roe.GetVector();
  tInit = tStart;
}
MatrixX RoeGeometricMappingDynamics::ComputeMatrix(real t) {
  real M = n * (t - tInit) + M0;
  real f = MeanToTrueAnomaly(M, e);
  real u = f + w;
  real sin_u = sin(u);
  real cos_u = cos(u);
  real cot_i = 1 / tan(i);
  real k = 1.0 + ex * cos(u) + ey * sin(u);
  real kp = -ex * sin(u) + ey * cos(u);
  real eta = sqrt(1.0 - e * e);

  MatrixX A = MatrixX::Zero(6, 6);
  A.block(0, 0, 3, 3) = a * eta * eta * MatrixXd::Identity(3, 3);
  A.block(3, 3, 3, 3) = a * n / eta * MatrixXd::Identity(3, 3);

  MatrixX B = MatrixX::Zero(6, 6);
  B(0, 0) = 1 / k + 3 / 2 * kp * n / pow(eta, 3) * t;
  B(0, 1) = -kp / pow(eta, 3);
  B(0, 2) = 1 / pow(eta, 3) * (ex * (k - 1) / (1 + eta) - cos_u);
  B(0, 3) = 1 / pow(eta, 3) * (ey * (k - 1) / (1 + eta) - sin_u);
  B(0, 5) = kp / pow(eta, 3) * cot_i;

  B(1, 0) = -3 / 2 * k * n / pow(eta, 3) * t;
  B(1, 1) = k / pow(eta, 3);
  B(1, 2) = 1 / pow(eta, 2) *
            ((1 + 1 / k) * sin_u + ey / k + k / eta * (ey / (1 + eta)));
  B(1, 3) = -1 / pow(eta, 2) *
            ((1 + 1 / k) * cos_u + ex / k + k / eta * (ex / (1 + eta)));
  B(1, 5) = (1 / k - k / pow(eta, 3)) * cot_i;

  B(2, 4) = 1 / k * sin_u;
  B(2, 5) = -1 / k * cos_u;

  B(3, 0) = kp / 2 + 3 / 2 * k * k * (1 - k) * n / pow(eta, 3) * t;
  B(3, 1) = k * k / pow(eta, 3) * (k - 1);
  B(3, 2) = k * k / pow(eta, 3) * (eta * sin_u + ey * (k - 1) / (1 + eta));
  B(3, 3) = -k * k / pow(eta, 3) * (eta * cos_u + ex * (k - 1) / (1 + eta));
  B(3, 5) = -k * k / pow(eta, 3) * (k - 1) * cot_i;

  B(4, 0) = -3 / 2 * k * (1 + k * kp * n / pow(eta, 3) * t);
  B(4, 1) = k * k / pow(eta, 3) * kp;
  B(4, 2) = (1 + k * k / pow(eta, 3)) * cos_u +
            ex * k / pow(eta, 2) * (1 + k / eta * (1 - k) / (1 + eta));
  B(4, 3) = (1 + k * k / pow(eta, 3)) * sin_u +
            ey * k / pow(eta, 2) * (1 + k / eta * (1 - k) / (1 + eta));
  B(4, 5) = -(1 + k * k / pow(eta, 3)) * kp * cot_i;

  B(5, 4) = cos_u + ex;
  B(5, 5) = sin_u + ey;

  MatrixX Phi = A * B;
  std::cout << A << std::endl << std::endl;
  std::cout << B << std::endl << std::endl;
  std::cout << Phi << std::endl;
  return Phi;
}

}  // namespace lupnt