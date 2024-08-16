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

#include "lupnt/core/constants.h"
#include "lupnt/core/progress_bar.h"
#include "lupnt/dynamics/dynamics.h"

namespace lupnt {

  // MatX IDynamics::Propagate(const MatX &x0, Real t0, Real tf) {
  //   MatX xf = MatX::Zero(x0.rows(), x0.cols());
  //   for (int i = 0; i < x0.cols(); i++) {
  //     Vec6 x0_col = x0.col(i);
  //     Vec6 xf_col = Propagate(x0_col, t0, tf);
  //     xf.col(i) = xf_col;
  //   }
  //   return xf;
  // }

  // MatX IDynamics::Propagate(const Vec6 &x0, Real t0, const Vec6 &tf, bool progress) {
  //   MatX xf = MatX::Zero(x0.rows(), tf.size());
  //   ProgressBar pbar(tf.size());
  //   for (int i = 0; i < tf.size(); i++) {
  //     Vec6 x0_col = x0;
  //     Real tf_i = tf(i);
  //     Vec6 xf_col = Propagate(x0_col, t0, tf_i);
  //     xf.col(i) = xf_col;
  //     if (progress) pbar.Update(i);
  //   }
  //   if (progress) pbar.Finish();
  //   return xf;
  // }

  NumericalOrbitDynamics::NumericalOrbitDynamics(ODE odefunc, IntegratorType integrator)
      : odefunc_(odefunc), propagator_(integrator) {}

  void NumericalOrbitDynamics::SetTimeStep(Real dt) { dt_ = dt.val(); };

  Vec6 NumericalOrbitDynamics::Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm) {
    Real dt_prop = (abs(dt_) < EPS) ? (tf - t0) / 10 : dt_;
    if (stm == nullptr) {
      Vec6 xf = propagator_.Propagate(odefunc_, t0, tf, x0, dt_prop);
      return xf;
    } else {
      MatXd stm_X(6, 6);
      Vec6 xf = propagator_.Propagate(odefunc_, t0, tf, x0, dt_prop, &stm_X);
      *stm = stm_X;
      return xf;
    }
  }

  MatX6 NumericalOrbitDynamics::Propagate(const Vec6 &x0, Real t0, const VecX &tf, bool progress) {
    MatX6 xf = MatX6::Zero(x0.rows(), tf.size());
    ProgressBar pbar(tf.size());
    xf.row(0) = Propagate(x0, t0, tf(0));
    for (int i = 1; i < tf.size(); i++) {
      Real t0_i = tf(i - 1);
      Real tf_i = tf(i);
      Vec6 xf_i = Propagate(x0, t0_i, tf_i);
      xf.row(i) = xf_i;
      if (progress) pbar.Update(i);
    }
    if (progress) pbar.Finish();
    return xf;
  }

  // ****************************************************************************
  // CartesianTwoBodyDynamics
  // ****************************************************************************

  CartesianTwoBodyDynamics::CartesianTwoBodyDynamics(Real GM, IntegratorType integ)
      : NumericalOrbitDynamics(std::bind(&CartesianTwoBodyDynamics::ComputeRates, this,
                                         std::placeholders::_1, std::placeholders::_2),
                               integ),
        GM_(GM) {};

  Vec6 CartesianTwoBodyDynamics::ComputeRates(Real t, const Vec6 &x) const {
    (void)t;
    Vec6 rv_dot;
    Vec3 r = x.head(3);
    Vec3 v = x.tail(3);
    Real r_norm = r.norm();

    rv_dot.head(3) = v;
    rv_dot.tail(3) = -GM_ * r / pow(r_norm, 3);

    return rv_dot;
  }

  // ****************************************************************************
  // J2CartTwoBodyDynamics
  // ****************************************************************************

  J2CartTwoBodyDynamics::J2CartTwoBodyDynamics(Real GM, Real J2, Real R_body, IntegratorType integ)
      : NumericalOrbitDynamics(std::bind(&J2CartTwoBodyDynamics::ComputeRates, this,
                                         std::placeholders::_1, std::placeholders::_2),
                               integ),
        GM_(GM),
        J2_(J2),
        R_body_(R_body) {};

  Vec6 J2CartTwoBodyDynamics::ComputeRates(Real t, const Vec6 &x) const {
    (void)t;
    Vec6 rv_dot(6);
    Vec3 r = x.head(3);
    Vec3 v = x.tail(3);
    Real r_norm = r.norm();

    rv_dot.head(3) = v;
    rv_dot.tail(3) = -GM_ * r / pow(r_norm, 3);

    Vec3 a_J2;
    Real aux1 = -3.0 / 2.0 * GM_ * J2_ * pow(R_body_, 2.0) / pow(r_norm, 4.0);
    Real aux2 = 5.0 * pow(r(2) / R_body_, 2.0);
    a_J2(0) = aux1 * (1.0 - aux2) * r(0);
    a_J2(1) = aux1 * (1.0 - aux2) * r(1);
    a_J2(2) = aux1 * (3.0 - aux2) * r(2);

    rv_dot.tail(3) += a_J2;

    return rv_dot;
  }

  // ****************************************************************************
  // J2KeplerianDynamics
  // ****************************************************************************

  J2KeplerianDynamics::J2KeplerianDynamics(Real GM, Real J2, Real R_body, IntegratorType integ)
      : NumericalOrbitDynamics(std::bind(&J2KeplerianDynamics::ComputeRates, this,
                                         std::placeholders::_1, std::placeholders::_2),
                               integ),
        GM_(GM),
        J2_(J2),
        R_body_(R_body) {};

  Vec6 J2KeplerianDynamics::ComputeRates(Real t, const Vec6 &x) const {
    (void)t;
    Real p = x(0) * (1.0 - x(1) * x(1));
    Real n = sqrt(GM_ / pow(x(0), 3.0));
    Real eta = sqrt(1.0 - x(1) * x(1));

    Vec6 coe_dot(6);
    coe_dot(0) = 0.0;
    coe_dot(1) = 0.0;
    coe_dot(2) = 0.0;
    coe_dot(3) = -3.0 / 2.0 * J2_ * pow(R_body_ / p, 2.0) * n * cos(x(2));
    coe_dot(4) = 3.0 / 4.0 * J2_ * pow(R_body_ / p, 2.0) * n * (5.0 * pow(cos(x(2)), 2.0) - 1.0);
    coe_dot(5)
        = n + 3.0 / 4.0 * J2_ * pow(R_body_ / p, 2.0) * n * eta * (3.0 * pow(cos(x(2)), 2.0) - 1.0);
    return coe_dot;
  }

  // ****************************************************************************
  // MoonMeanDynamics
  // ****************************************************************************

  MoonMeanDynamics::MoonMeanDynamics(IntegratorType integ)
      : NumericalOrbitDynamics(std::bind(&MoonMeanDynamics::ComputeRates, this,
                                         std::placeholders::_1, std::placeholders::_2),
                               integ) {};

  Vec6 MoonMeanDynamics::ComputeRates(Real t, const Vec6 &x) const {
    (void)t;
    Real a = x(0);
    Real e = x(1);
    Real i = x(2);
    Real w = x(4);

    Vec6 coe_dot(6);
    coe_dot(0) = 0;
    coe_dot(1) = (15 * k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0)) / (8 * sqrt(GM_MOON)) * e
                 * sqrt(1 - pow(e, 2)) * pow(sin(i), 2) * sin(2 * w);
    coe_dot(2) = -(15 * k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0)) / (16 * sqrt(GM_MOON)) * pow(e, 2)
                 / sqrt(1 - pow(e, 2)) * sin(2 * i) * sin(2 * w);
    coe_dot(3) = -(3 * J2_ * sqrt(GM_MOON) * pow(R_MOON, 2))
                     / (2 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 2)) * cos(i)
                 + (3 * k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0))
                       / (8 * sqrt(GM_MOON) * sqrt(1 - pow(e, 2)))
                       * (5 * pow(e, 2) * cos(2 * w) - 3 * pow(e, 2) - 2) * cos(i);
    coe_dot(4) = (3 * J2_ * sqrt(GM_MOON) * pow(R_MOON, 2))
                     / (4 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 2)) * (5 * pow(cos(i), 2) - 1)
                 + (3 * k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0))
                       / (8 * sqrt(GM_MOON) * sqrt(1 - pow(e, 2)))
                       * ((5 * pow(cos(i), 2) - 1 + pow(e, 2))
                          + 5 * (1 - pow(e, 2) - pow(cos(i), 2)) * cos(2 * w));
    coe_dot(5) = sqrt(GM_MOON) / pow(a, 3.0 / 2.0)
                 + (3 * J2_ * sqrt(GM_MOON) * pow(R_MOON, 2))
                       / (4 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 3.0 / 2.0))
                       * (3 * pow(cos(i), 2) - 1)
                 - (k_ * pow(n3_, 2) * pow(a, 3.0 / 2.0)) / (8 * sqrt(GM_MOON))
                       * ((3 * pow(e, 2) + 7) * (3 * pow(cos(i), 2) - 1)
                          + 15 * (1 + pow(e, 2)) * pow(sin(i), 2) * cos(2 * w));
    return coe_dot;
  }
};  // namespace lupnt
