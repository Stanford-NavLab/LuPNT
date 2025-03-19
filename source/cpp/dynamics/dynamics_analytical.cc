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
#include <lupnt/core/constants.h>
#include <lupnt/core/progress_bar.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/orbit_state.h>

#include <Eigen/QR>

namespace lupnt {

  MatX6 IAnalyticalOrbitDynamics::Propagate(const Vec6 &x0, Real t0, const VecX &tf,
                                            bool progress) {
    MatX6 xf = MatX6::Zero(tf.size(), 6);
    ProgressBar pbar(tf.size());
    for (int i = 0; i < tf.size(); i++) {
      Real tf_i = tf(i);
      Vec6 xf_i = Propagate(x0, t0, tf_i);
      xf.row(i) = xf_i;
      if (progress) pbar.Update(i);
    }
    if (progress) pbar.Finish();
    return xf;
  }

  // ****************************************************************************
  // KeplerianDynamics
  // ****************************************************************************

  KeplerianDynamics::KeplerianDynamics(Real GM) : GM_(GM) {};

  OrbitState KeplerianDynamics::PropagateState(const OrbitState &state, Real t0, Real tf,
                                               Mat6d *stm) {
    switch (state.GetOrbitStateRepres()) {
      case OrbitStateRepres::CLASSICAL_OE: {
        Vec6 xf = PropagateClassicalOE(state.GetVec(), t0, tf, stm);
        return ClassicalOE(xf, state.GetFrame());
      }
      case OrbitStateRepres::QUASI_NONSINGULAR_OE: {
        Vec6 xf = PropagateQuasiNonsingOE(state.GetVec(), t0, tf, stm);
        return QuasiNonsingOE(xf, state.GetFrame());
      }
      case OrbitStateRepres::EQUINOCTIAL_OE: {
        Vec6 xf = PropagateEquinoctialOE(state.GetVec(), t0, tf, stm);
        return EquinoctialOE(xf, state.GetFrame());
      }
      default: break;
    }
    throw std::runtime_error("OrbitState type not supported");
  }

  Vec6 KeplerianDynamics::Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm) {
    if (abs(tf - t0) < EPS) {
      if (stm != nullptr) *stm = Mat6d::Identity(6, 6);
      return x0;
    }
    return PropagateClassicalOE(x0, t0, tf, stm);
  }

  // ClassicalOE
  Vec6 KeplerianDynamics::PropagateClassicalOE(const Vec6 &coe, Real t0, Real tf, Mat6d *stm) {
    Real dt = tf - t0;
    Real n = sqrt(GM_ / pow(coe[0], 3));
    Real M = Wrap2Pi(coe[5] + n * dt);

    Vec6 coe_new = coe;
    coe_new[5] = M;

    if (stm != nullptr) {
      *stm = Mat6d::Identity(6, 6);
      (*stm)(5, 0) = -3. / 2. * (n / coe[0] * dt).val();
    }
    return coe_new;
  }

  // QuasiNonsingOE
  Vec6 KeplerianDynamics::PropagateQuasiNonsingOE(const Vec6 &qnsoe, Real t0, Real tf, Mat6d *stm) {
    Real dt = tf - t0;
    Real n = sqrt(GM_ / pow(qnsoe[0], 3));
    Real u = qnsoe[1] + n * dt;

    Vec6 qnsoe_new = qnsoe;
    qnsoe_new[1] = u;

    if (stm != nullptr) throw std::runtime_error("Not implemented");
    return qnsoe_new;
  }

  // EquinoctialOE
  Vec6 KeplerianDynamics::PropagateEquinoctialOE(const Vec6 &eqoe, Real t0, Real tf, Mat6d *stm) {
    Real dt = tf - t0;
    Real n = sqrt(GM_ / pow(eqoe[0], 3));
    Real lon = eqoe[5] + n * dt;

    Vec6 eqoe_new = eqoe;
    eqoe_new[5] = lon;

    if (stm != nullptr) throw std::runtime_error("Not implemented");
    return eqoe_new;
  }

  /* ****************************************************************************
    ClohessyWiltshireDynamics
    ************************************************************************** */

  ClohessyWiltshireDynamics::ClohessyWiltshireDynamics(Real a, Real n) : a_(a), n_(n) {};

  Vec6 ClohessyWiltshireDynamics::Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm) {
    if (abs(tf - t0) < EPS) {
      if (stm != nullptr) *stm = Mat6d::Identity(6, 6);
      return x0;
    }
    if (t0 != t0_) {
      MatX Phi = ComputeMat(t0);
      K_ = Phi.colPivHouseholderQr().solve(x0);
    }
    Mat6 Phi = ComputeMat(tf - t0);
    Vec6 xf = Phi * K_;
    if (stm != nullptr) *stm = Phi.cast<double>();
    return xf;
  }

  OrbitState ClohessyWiltshireDynamics::PropagateState(const OrbitState &state, Real t0, Real tf,
                                                       Mat6d *stm) {
    assert(state.GetOrbitStateRepres() == OrbitStateRepres::CARTESIAN
           && "OrbitState type not supported");
    Vec6 xf = Propagate(state.GetVec(), t0, tf, stm);
    return CartesianOrbitState(xf, state.GetFrame());
  }

  Mat6 ClohessyWiltshireDynamics::ComputeMat(Real t) {
    Real sin_nt = sin(n_ * t);
    Real cos_nt = cos(n_ * t);

    Mat6 A = Mat6::Zero();
    A.block(0, 0, 3, 3) = a_ * Vec6d::Identity(3, 3);
    A.block(3, 3, 3, 3) = a_ * n_ * Vec6d::Identity(3, 3);

    Mat6 B = Mat6::Zero();
    B(0, 0) = 1.;
    B(0, 1) = sin_nt;
    B(0, 2) = cos_nt;

    B(1, 0) = -3. / 2. * n_ * t;
    B(1, 1) = 2. * cos_nt;
    B(1, 2) = -2. * sin_nt;
    B(1, 3) = 1.;

    B(2, 4) = sin_nt;
    B(2, 5) = cos_nt;

    B(3, 1) = cos_nt;
    B(3, 2) = -sin_nt;

    B(4, 0) = -3. / 2.;
    B(4, 1) = -2. * sin_nt;
    B(4, 2) = -2. * cos_nt;

    B(5, 4) = cos_nt;
    B(5, 5) = -sin_nt;

    Mat6 Phi = A * B;
    return Phi;
  }

  // ****************************************************************************
  // YamanakaAnkersenDynamics
  // ****************************************************************************

  YamanakaAnkersenDynamics::YamanakaAnkersenDynamics(const ClassicalOE &coe_c,
                                                     const CartesianOrbitState &rv_rtn, Real GM_) {
    a_ = coe_c.a().val();
    n_ = sqrt(GM_ / pow(a_, 3.));
    e_ = coe_c.e().val();
    M0_ = coe_c.M().val();
    rv_rtn_ = rv_rtn.GetVec();
  }

  Vec6 YamanakaAnkersenDynamics::Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm) {
    if (abs(tf - t0) < EPS) {
      if (stm != nullptr) *stm = Mat6d::Identity(6, 6);
      return x0;
    }
    throw std::runtime_error("Not implemented");

    if (t0 != t0_) {
      MatX Phi = ComputeMat(t0);
      K_ = ComputeInverseMat(t0) * rv_rtn_;
      // K_ = Phi.colPivHouseholderQr().solve(x0);
    }
    Mat6 Phi = ComputeMat(tf - t0);
    Vec6 xf = Phi * K_;
    if (stm != nullptr) *stm = Phi.cast<double>();
    return xf;
  }

  OrbitState YamanakaAnkersenDynamics::PropagateState(const OrbitState &state, Real t0, Real tf,
                                                      Mat6d *stm) {
    assert(state.GetOrbitStateRepres() == OrbitStateRepres::CARTESIAN
           && "OrbitState type not supported");
    Vec6 xf = Propagate(state.GetVec(), t0, tf, stm);
    return CartesianOrbitState(xf, state.GetFrame());
  }

  MatX YamanakaAnkersenDynamics::ComputeMat(Real t) {
    Real M = n_ * (t - t0_) + M0_;
    Real f = Mean2TrueAnomaly(M, e_);
    Real sin_f = sin(f);
    Real cos_f = cos(f);
    Real k = 1. + e_ * cos(f);
    Real kp = -e_ * sin(f);
    Real eta = sqrt(1. - e_ * e_);
    Real tau = n_ * t / pow(eta, 3.);

    MatX A = MatX::Zero(6, 6);
    A.block(0, 0, 3, 3) = a_ * eta * eta * Vec6d::Identity(3, 3);
    A.block(3, 3, 3, 3) = a_ * n_ / eta * Vec6d::Identity(3, 3);

    MatX B = MatX::Zero(6, 6);
    B(0, 0) = 1. / k + 3. / 2. * kp * tau;
    B(0, 1) = sin_f;
    B(0, 2) = cos_f;

    B(1, 0) = -3. / 2. * k * tau;
    B(1, 1) = (1. + 1. / k) * cos_f;
    B(1, 2) = -(1. + 1. / k) * sin_f;
    B(1, 3) = 1. / k;

    B(2, 4) = 1. / k * sin_f;
    B(2, 5) = 1. / k * cos_f;

    B(3, 0) = kp / 2. - 3. / 2. * k * k * (k - 1.) * tau;
    B(3, 1) = k * k * cos_f;
    B(3, 2) = -k * k * sin_f;

    B(4, 0) = -3. / 2. * (k + k * k * kp * tau);
    B(4, 1) = -(k * k + 1) * sin_f;
    B(4, 2) = -e_ - (k * k + 1) * cos_f;
    B(4, 3) = -kp;

    B(5, 4) = e_ + cos_f;
    B(5, 5) = -sin_f;

    MatX Phi = A * B;
    return Phi;
  }
  MatX YamanakaAnkersenDynamics::ComputeInverseMat(Real t) {
    Real M = n_ * (t - t0_) + M0_;
    Real f = Mean2TrueAnomaly(M, e_);
    Real sin_f = sin(f);
    Real cos_f = cos(f);
    Real k = 1. + e_ * cos(f);
    Real kp = -e_ * sin(f);
    Real eta = sqrt(1. - e_ * e_);
    Real tau = n_ * t / pow(eta, 3.);

    MatX A = MatX::Zero(6, 6);
    A.block(0, 0, 3, 3) = 1. / a_ / (eta * eta) * Vec6d::Identity(3, 3);
    A.block(3, 3, 3, 3) = eta / a_ / n_ * Vec6d::Identity(3, 3);

    MatX B = MatX::Zero(6, 6);
    B(0, 0) = 2. * (k * k) * (k + 1) / (eta * eta);
    B(0, 1) = 2. * (k * k) * kp / (eta * eta);
    B(0, 3) = -2. * kp / (eta * eta);
    B(0, 4) = 2. * k / (eta * eta);

    B(1, 0) = (1. - pow(k + 1, 2) / (eta * eta)) * sin_f
              + 3. * e_ * (k * k) * (k + 1) * tau / (eta * eta);
    B(1, 1) = -(k + 1) * kp * sin_f / (eta * eta) + 3. * e_ * (k * k) * kp * tau / (eta * eta);
    B(1, 3) = (1. / (eta * eta)) * (cos_f - 2. * e_ / k) - 3. * e_ * kp * tau / (eta * eta);
    B(1, 4) = -(1. / (eta * eta)) * (1. + 1. / k) * sin_f + 3. * e_ * k * tau / (eta * eta);

    B(2, 0) = -(k / (eta * eta)) * (2. * e_ + (k + 2) * cos_f);
    B(2, 1) = -(1. / (eta * eta)) * (e_ + (k + 1) * cos_f) * kp;
    B(2, 3) = -(1. / (eta * eta)) * sin_f;
    B(2, 4) = -(1. / (eta * eta)) * (e_ / k + (1. + 1. / k) * cos_f);

    B(3, 0) = (pow(k + 1, 2) / (eta * eta)) * kp + 3. * (k * k) * (k + 1) * tau / (eta * eta);
    B(3, 1) = (k / (eta * eta)) * (2. + k - (k * k)) - 1. + 3. * (k * k) * kp * tau / (eta * eta);
    B(3, 3) = (1. / (eta * eta)) * (k - 1. - 2. / k) - 3. * kp * tau / (eta * eta);
    B(3, 4) = (1. / (eta * eta)) * (1. + 1. / k) * kp + 3. * k * tau / (eta * eta);

    B(4, 2) = sin_f;
    B(4, 5) = (1. / k) * cos_f;

    B(5, 2) = e_ + cos_f;
    B(5, 5) = -(1. / k) * sin_f;

    MatX Phi = B * A;
    return Phi;
  }

  // ****************************************************************************
  // RoeGeometricMappingDynamics
  // ****************************************************************************

  RoeGeometricMappingDynamics::RoeGeometricMappingDynamics(const ClassicalOE coe_c,
                                                           const QuasiNonsingROE &roe, Real GM) {
    a_ = coe_c.a();
    e_ = coe_c.e();
    i_ = coe_c.i();
    w_ = coe_c.w();
    M0_ = coe_c.M();

    ex_ = coe_c.e() * cos(coe_c.w());
    ey_ = coe_c.e() * sin(coe_c.w());

    n_ = sqrt(GM / pow(coe_c.a(), 3));
    K_ = roe.GetVec();
  }

  Vec6 RoeGeometricMappingDynamics::Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm) {
    if (abs(tf - t0) < EPS) {
      if (stm != nullptr) *stm = Mat6d::Identity(6, 6);
      return x0;
    }
    throw std::runtime_error("Not implemented");

    if (t0 != t0_) {
      throw std::runtime_error("Not implemented");
    }

    MatX Phi = ComputeMat(tf - t0_);
    Vec6 xf = Phi * K_;
    if (stm != nullptr) *stm = Phi.cast<double>();
    return xf;
  }

  OrbitState RoeGeometricMappingDynamics::PropagateState(const OrbitState &state, Real t0, Real tf,
                                                         Mat6d *stm) {
    assert(state.GetOrbitStateRepres() == OrbitStateRepres::CARTESIAN
           && "OrbitState type not supported");
    Vec6 xf = Propagate(state.GetVec(), t0, tf, stm);
    return CartesianOrbitState(xf, state.GetFrame());
  }

  MatX RoeGeometricMappingDynamics::ComputeMat(Real t) {
    Real M = n_ * (t - t0_) + M0_;
    Real f = Mean2TrueAnomaly(M, e_);
    Real u = f + w_;
    Real sin_u = sin(u);
    Real cos_u = cos(u);
    Real cot_i = 1. / tan(i_);
    Real k = 1. + ex_ * cos(u) + ey_ * sin(u);
    Real kp = -ex_ * sin(u) + ey_ * cos(u);
    Real eta = sqrt(1. - e_ * e_);

    MatX A = MatX::Zero(6, 6);
    A.block(0, 0, 3, 3) = a_ * eta * eta * Vec6d::Identity(3, 3);
    A.block(3, 3, 3, 3) = a_ * n_ / eta * Vec6d::Identity(3, 3);

    MatX B = MatX::Zero(6, 6);
    B(0, 0) = 1. / k + 3. / 2. * kp * n_ / pow(eta, 3) * t;
    B(0, 1) = -kp / pow(eta, 3);
    B(0, 2) = 1. / pow(eta, 3) * (ex_ * (k - 1.) / (1. + eta) - cos_u);
    B(0, 3) = 1. / pow(eta, 3) * (ey_ * (k - 1.) / (1. + eta) - sin_u);
    B(0, 5) = kp / pow(eta, 3) * cot_i;

    B(1, 0) = -3 / 2. * k * n_ / pow(eta, 3) * t;
    B(1, 1) = k / pow(eta, 3);
    B(1, 2) = 1. / pow(eta, 2) * ((1. + 1. / k) * sin_u + ey_ / k + k / eta * (ey_ / (1. + eta)));
    B(1, 3) = -1. / pow(eta, 2) * ((1. + 1. / k) * cos_u + ex_ / k + k / eta * (ex_ / (1. + eta)));
    B(1, 5) = (1. / k - k / pow(eta, 3)) * cot_i;

    B(2, 4) = 1. / k * sin_u;
    B(2, 5) = -1. / k * cos_u;

    B(3, 0) = kp / 2. + 3. / 2. * k * k * (1. - k) * n_ / pow(eta, 3) * t;
    B(3, 1) = k * k / pow(eta, 3) * (k - 1.);
    B(3, 2) = k * k / pow(eta, 3) * (eta * sin_u + ey_ * (k - 1.) / (1. + eta));
    B(3, 3) = -k * k / pow(eta, 3) * (eta * cos_u + ex_ * (k - 1.) / (1. + eta));
    B(3, 5) = -k * k / pow(eta, 3) * (k - 1.) * cot_i;

    B(4, 0) = -3 / 2. * k * (1. + k * kp * n_ / pow(eta, 3) * t);
    B(4, 1) = k * k / pow(eta, 3) * kp;
    B(4, 2) = (1. + k * k / pow(eta, 3)) * cos_u
              + ex_ * k / pow(eta, 2) * (1. + k / eta * (1. - k) / (1. + eta));
    B(4, 3) = (1. + k * k / pow(eta, 3)) * sin_u
              + ey_ * k / pow(eta, 2) * (1. + k / eta * (1. - k) / (1. + eta));
    B(4, 5) = -(1. + k * k / pow(eta, 3)) * kp * cot_i;

    B(5, 4) = cos_u + ex_;
    B(5, 5) = sin_u + ey_;

    MatX Phi = A * B;
    std::cout << A << std::endl << std::endl;
    std::cout << B << std::endl << std::endl;
    std::cout << Phi << std::endl;
    return Phi;
  }

}  // namespace lupnt
