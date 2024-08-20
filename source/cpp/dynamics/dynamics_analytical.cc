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
#include <Eigen/QR>

#include "lupnt/core/constants.h"
#include "lupnt/dynamics/dynamics.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/orbit_state.h"

namespace lupnt {

  // ****************************************************************************
  // KeplerianDynamics
  // ****************************************************************************

  KeplerianDynamics::KeplerianDynamics(double GM) : mu_(GM) {};

  // ClassicalOE
  void KeplerianDynamics::Propagate(ClassicalOE &state, Real dt) {
    Real a = state.a();
    Real n = sqrt(mu_ / pow(a, 3));
    state.Set_M(Wrap2Pi(state.M() + n * dt));
  }

  Vec6 KeplerianDynamics::PropagateClassicalOE(Vec6 coe, Real dt, double GM) {
    Real a = coe[0];
    Real n = sqrt(GM / pow(a, 3));
    coe[5] = Wrap2Pi(coe[5] + n * dt);
    return coe;
  }

  void KeplerianDynamics::PropagateWithStm(ClassicalOE &state, Real dt, Mat6d &stm) {
    Real a = state.a();
    Real n = sqrt(mu_ / pow(a, 3));
    state.Set_M(Wrap2Pi(state.M() + n * dt));
    stm = Mat6d::Identity(6, 6);
    stm(5, 0) = -3.0 / 2.0 * (n / a * dt).val();
  }

  // QuasiNonsingOE
  void KeplerianDynamics::Propagate(QuasiNonsingOE &state, Real dt) {
    state.Set_u(state.u() + sqrt(mu_ / pow(state.a(), 3)) * dt);
  }
  void KeplerianDynamics::PropagateWithStm(QuasiNonsingOE &state, Real dt, Mat6d &stm) {
    throw std::runtime_error("Not implemented");
  }

  // EquinoctialOE
  void KeplerianDynamics::Propagate(EquinoctialOE &state, Real dt) {
    state.Set_lon(state.lon() + sqrt(mu_ / pow(state.a(), 3)) * dt);
  }
  void KeplerianDynamics::PropagateWithStm(EquinoctialOE &state, Real dt, Mat6d &stm) {
    throw std::runtime_error("Not implemented");
  }

  /* ****************************************************************************
    ClohessyWiltshireDynamics
    ************************************************************************** */

  ClohessyWiltshireDynamics::ClohessyWiltshireDynamics(Real a_in, Real n_in) : a(a_in), n(n_in) {};

  void ClohessyWiltshireDynamics::Propagate(OrbitState &state, Real tEnd) {
    if (state.GetOrbitStateRepres() != OrbitStateRepres::CARTESIAN)
      throw std::runtime_error("OrbitState type not supported");

    VecX xEnd = ComputeMat(tEnd) * K;
    state.SetVec(xEnd);
  }

  void ClohessyWiltshireDynamics::Initialize(CartesianOrbitState &state, Real tStart) {
    tInit = tStart;
    MatX Phi = ComputeMat(tStart);
    K = Phi.colPivHouseholderQr().solve(state.GetVec());
  }

  MatX ClohessyWiltshireDynamics::ComputeMat(Real t) {
    Real sin_nt = sin(n * t);
    Real cos_nt = cos(n * t);

    MatX A = MatX::Zero(6, 6);
    A.block(0, 0, 3, 3) = a * VecXd::Identity(3, 3);
    A.block(3, 3, 3, 3) = a * n * VecXd::Identity(3, 3);

    MatX B = MatX::Zero(6, 6);
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

    MatX Phi = A * B;
    return Phi;
  }

  // ****************************************************************************
  // YamanakaAnkersenDynamics
  // ****************************************************************************

  YamanakaAnkersenDynamics::YamanakaAnkersenDynamics() : a(0.0), n(0.0), e(0.0), M0(0.0) {};
  void YamanakaAnkersenDynamics::Propagate(CartesianOrbitState &state, Real tEnd) {
    if (state.GetOrbitStateRepres() == OrbitStateRepres::CARTESIAN) {
      VecX xEnd = ComputeMat(tEnd) * K;
      state.SetVec(xEnd);
    } else {
      throw std::runtime_error("OrbitState type not supported");
    }
  }
  void YamanakaAnkersenDynamics::Initialize(ClassicalOE &coe_c, CartesianOrbitState &rv_rtn,
                                            Real tStart, double GM) {
    a = coe_c.a().val();
    n = sqrt(GM / pow(a, 3.0));
    e = coe_c.e().val();
    M0 = coe_c.M().val();
    tInit = tStart;

    MatX Phi = ComputeMat(tStart);
    K = ComputeInverseMat(tStart) * rv_rtn.GetVec();
    // K = Phi.colPivHouseholderQr().solve(state.GetVec());
  }
  MatX YamanakaAnkersenDynamics::ComputeMat(Real t) {
    Real M = n * (t - tInit) + M0;
    Real f = Mean2TrueAnomaly(M, e);
    Real sin_f = sin(f);
    Real cos_f = cos(f);
    Real k = 1.0 + e * cos(f);
    Real kp = -e * sin(f);
    Real eta = sqrt(1.0 - e * e);
    Real tau = n * t / pow(eta, 3.0);

    MatX A = MatX::Zero(6, 6);
    A.block(0, 0, 3, 3) = a * eta * eta * VecXd::Identity(3, 3);
    A.block(3, 3, 3, 3) = a * n / eta * VecXd::Identity(3, 3);

    MatX B = MatX::Zero(6, 6);
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

    MatX Phi = A * B;
    return Phi;
  }
  MatX YamanakaAnkersenDynamics::ComputeInverseMat(Real t) {
    Real M = n * (t - tInit) + M0;
    Real f = Mean2TrueAnomaly(M, e);
    Real sin_f = sin(f);
    Real cos_f = cos(f);
    Real k = 1.0 + e * cos(f);
    Real kp = -e * sin(f);
    Real eta = sqrt(1.0 - e * e);
    Real tau = n * t / pow(eta, 3.0);

    MatX A = MatX::Zero(6, 6);
    A.block(0, 0, 3, 3) = 1.0 / a / (eta * eta) * VecXd::Identity(3, 3);
    A.block(3, 3, 3, 3) = eta / a / n * VecXd::Identity(3, 3);

    MatX B = MatX::Zero(6, 6);
    B(0, 0) = 2.0 * (k * k) * (k + 1) / (eta * eta);
    B(0, 1) = 2.0 * (k * k) * kp / (eta * eta);
    B(0, 3) = -2.0 * kp / (eta * eta);
    B(0, 4) = 2.0 * k / (eta * eta);

    B(1, 0) = (1.0 - pow(k + 1, 2) / (eta * eta)) * sin_f
              + 3.0 * e * (k * k) * (k + 1) * tau / (eta * eta);
    B(1, 1) = -(k + 1) * kp * sin_f / (eta * eta) + 3.0 * e * (k * k) * kp * tau / (eta * eta);
    B(1, 3) = (1.0 / (eta * eta)) * (cos_f - 2.0 * e / k) - 3.0 * e * kp * tau / (eta * eta);
    B(1, 4) = -(1.0 / (eta * eta)) * (1.0 + 1.0 / k) * sin_f + 3.0 * e * k * tau / (eta * eta);

    B(2, 0) = -(k / (eta * eta)) * (2.0 * e + (k + 2) * cos_f);
    B(2, 1) = -(1.0 / (eta * eta)) * (e + (k + 1) * cos_f) * kp;
    B(2, 3) = -(1.0 / (eta * eta)) * sin_f;
    B(2, 4) = -(1.0 / (eta * eta)) * (e / k + (1.0 + 1.0 / k) * cos_f);

    B(3, 0) = (pow(k + 1, 2) / (eta * eta)) * kp + 3.0 * (k * k) * (k + 1) * tau / (eta * eta);
    B(3, 1)
        = (k / (eta * eta)) * (2.0 + k - (k * k)) - 1.0 + 3.0 * (k * k) * kp * tau / (eta * eta);
    B(3, 3) = (1.0 / (eta * eta)) * (k - 1.0 - 2.0 / k) - 3.0 * kp * tau / (eta * eta);
    B(3, 4) = (1.0 / (eta * eta)) * (1.0 + 1.0 / k) * kp + 3.0 * k * tau / (eta * eta);

    B(4, 2) = sin_f;
    B(4, 5) = (1.0 / k) * cos_f;

    B(5, 2) = e + cos_f;
    B(5, 5) = -(1.0 / k) * sin_f;

    MatX Phi = B * A;
    return Phi;
  }

  // ****************************************************************************
  // RoeGeometricMappingDynamics
  // ****************************************************************************

  RoeGeometricMappingDynamics::RoeGeometricMappingDynamics()
      : a(0.0), n(0.0), e(0.0), M0(0.0), ex(0.0), ey(0.0), tInit(0.0) {};
  void RoeGeometricMappingDynamics::Propagate(CartesianOrbitState &state, Real tEnd) {
    if (state.GetOrbitStateRepres() == OrbitStateRepres::CARTESIAN) {
      std::cout << "xStart: " << state.GetVec().transpose() << std::endl;
      VecX xEnd = ComputeMat(tEnd) * K;
      state.SetVec(xEnd);
      std::cout << "K: " << K.transpose() << std::endl;
      std::cout << "xEnd: " << xEnd.transpose() << std::endl;
    } else {
      throw std::runtime_error("OrbitState type not supported");
    }
  }
  void RoeGeometricMappingDynamics::Initialize(ClassicalOE coe_c, QuasiNonsingROE &roe, Real tStart,
                                               double GM) {
    a = coe_c.a().val();
    e = coe_c.e().val();
    i = coe_c.i().val();
    w = coe_c.w().val();
    M0 = coe_c.M().val();

    ex = (double)(coe_c.e() * cos(coe_c.w()));
    ey = (double)(coe_c.e() * sin(coe_c.w()));

    n = (double)sqrt(GM / pow(coe_c.a(), 3));
    K = roe.GetVec();
    tInit = tStart;
  }
  MatX RoeGeometricMappingDynamics::ComputeMat(Real t) {
    Real M = n * (t - tInit) + M0;
    Real f = Mean2TrueAnomaly(M, e);
    Real u = f + w;
    Real sin_u = sin(u);
    Real cos_u = cos(u);
    Real cot_i = 1 / tan(i);
    Real k = 1.0 + ex * cos(u) + ey * sin(u);
    Real kp = -ex * sin(u) + ey * cos(u);
    Real eta = sqrt(1.0 - e * e);

    MatX A = MatX::Zero(6, 6);
    A.block(0, 0, 3, 3) = a * eta * eta * VecXd::Identity(3, 3);
    A.block(3, 3, 3, 3) = a * n / eta * VecXd::Identity(3, 3);

    MatX B = MatX::Zero(6, 6);
    B(0, 0) = 1 / k + 3 / 2 * kp * n / pow(eta, 3) * t;
    B(0, 1) = -kp / pow(eta, 3);
    B(0, 2) = 1 / pow(eta, 3) * (ex * (k - 1) / (1 + eta) - cos_u);
    B(0, 3) = 1 / pow(eta, 3) * (ey * (k - 1) / (1 + eta) - sin_u);
    B(0, 5) = kp / pow(eta, 3) * cot_i;

    B(1, 0) = -3 / 2 * k * n / pow(eta, 3) * t;
    B(1, 1) = k / pow(eta, 3);
    B(1, 2) = 1 / pow(eta, 2) * ((1 + 1 / k) * sin_u + ey / k + k / eta * (ey / (1 + eta)));
    B(1, 3) = -1 / pow(eta, 2) * ((1 + 1 / k) * cos_u + ex / k + k / eta * (ex / (1 + eta)));
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
    B(4, 2) = (1 + k * k / pow(eta, 3)) * cos_u
              + ex * k / pow(eta, 2) * (1 + k / eta * (1 - k) / (1 + eta));
    B(4, 3) = (1 + k * k / pow(eta, 3)) * sin_u
              + ey * k / pow(eta, 2) * (1 + k / eta * (1 - k) / (1 + eta));
    B(4, 5) = -(1 + k * k / pow(eta, 3)) * kp * cot_i;

    B(5, 4) = cos_u + ex;
    B(5, 5) = sin_u + ey;

    MatX Phi = A * B;
    std::cout << A << std::endl << std::endl;
    std::cout << B << std::endl << std::endl;
    std::cout << Phi << std::endl;
    return Phi;
  }

}  // namespace lupnt
