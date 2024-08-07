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

  /* **************************************************************************
   *  Numerical Dynamics
   * **************************************************************************/

  // arbitrary state size
  void NumericalDynamics::PropagateX(VecX &x, Real t0, Real tf) {
    Real dt_prop = tf - t0;
    if (dt_ == 0.0) {
      dt_prop = (tf - t0) / 10;
    } else {
      dt_prop = dt_;
    }

    VecX x0 = x.cast<double>();  // cut the link with previous relations
    VecX xf = propagator_.Propagate(odefunc_, t0, tf, x0, dt_prop);
    x = xf;
  }

  void NumericalDynamics::PropagateWithStmX(VecX &x, Real t0, Real tf, MatXd &stm) {
    Real dt_prop = tf - t0;
    if (dt_ == 0.0) {
      dt_prop = (tf - t0) / 10;
    } else {
      dt_prop = dt_;
    }

    VecX x0 = x.cast<double>();  // cut the link with previous relations
    stm.resize(x.size(), x.size());
    VecX xf = propagator_.PropagateWithStm(odefunc_, t0, tf, x0, dt_prop, stm);
    x = xf;
  }

  // ****************************************************************************
  // NumericalOrbitDynamics
  // ****************************************************************************

  void NumericalOrbitDynamics::Propagate(OrbitState &state, Real t0, Real tf, Real dt) {
    assert(state.GetOrbitStateRepres() == state_representation_);
    Real dt_prop = (dt > 0.0) ? dt : dt_;
    Vec6 x0 = state.GetVec();
    Vec6 xf = propagator_.Propagate(odefunc_, t0, tf, x0, dt_prop);
    state.SetVec(xf);
  }

  void NumericalOrbitDynamics::Propagate(Vec6 &x, Real t0, Real tf, Real dt) {
    Real dt_prop = (dt > 0.0) ? dt : dt_;
    Vec6 x0 = x;
    Vec6 xf = propagator_.Propagate(odefunc_, t0, tf, x0, dt_prop);
    x = xf;
  }

  void NumericalOrbitDynamics::PropagateWithStm(OrbitState &state, Real t0, Real tf, Real dt,
                                                Mat6d &stm) {
    assert(state.GetOrbitStateRepres() == state_representation_);
    Vec6 x0 = state.GetVec();
    MatXd J(6, 6);
    VecX xf = propagator_.PropagateWithStm(odefunc_, t0, tf, x0, dt, J);
    state.SetVec(xf);
    stm = J;
  }

  void NumericalOrbitDynamics::PropagateWithStm(Vec6 &x, Real t0, Real tf, Real dt, Mat6d &stm) {
    Real dt_prop = (dt > 0.0) ? dt : dt_;
    Vec6 x0 = x;
    MatXd J(6, 6);
    VecX xf = propagator_.PropagateWithStm(odefunc_, t0, tf, x0, dt_prop, J);
    x = xf;
    stm = J;
  }

  MatX NumericalOrbitDynamics::Propagate(Vec6 x, Real t0, VecX tf, Real dt, bool progress) {
    Real dt_prop = (dt > 0.0) ? dt : dt_;
    Vec6 x0 = x;
    MatX xf(tf.size(), 6);

    ProgressBar pb(tf.size());
    for (int i = 0; i < tf.size(); i++) {
      if (i == 0)
        xf.row(i) = propagator_.Propagate(odefunc_, t0, tf(i), x0, dt_prop);
      else
        xf.row(i) = propagator_.Propagate(odefunc_, tf(i - 1), tf(i), xf.row(i - 1), dt_prop);
      if (progress) pb.Update();
    }
    return xf;
  }

  MatX NumericalOrbitDynamics::Propagate(OrbitState &state, Real t0, VecX &tf, Real dt,
                                         bool progress) {
    Vec6 x0 = state.GetVec();
    return Propagate(x0, t0, tf, progress);
  }

  void NumericalOrbitDynamics::PropagateWithStm(OrbitState &state, Real t0, Real tf, Mat6d &stm) {
    Real dt_prop = (dt_ > 0.0) ? dt_ : (tf - t0) / 10;
    PropagateWithStm(state, t0, tf, dt_prop, stm);
  }

  void NumericalOrbitDynamics::PropagateWithStm(Vec6 &x, Real t0, Real tf, Mat6d &stm) {
    Real dt_prop = (dt_ > 0.0) ? dt_ : (tf - t0) / 10;
    PropagateWithStm(x, t0, tf, dt_prop, stm);
  }
  // With returns --------------------------------------------------------
  Vec6 NumericalOrbitDynamics::PropagateR(Vec6 &x, Real t0, Real tf, Real dt) {
    Vec6 x0 = x;
    Vec6 xf = propagator_.Propagate(odefunc_, t0, tf, x0, dt);
    return xf;
  }

  Vec6 NumericalOrbitDynamics::PropagateWithStmR(Vec6 &x, Real t0, Real tf, Real dt, Mat6d &stm) {
    Vec6 x0 = x;
    MatXd J(6, 6);
    VecX xf = propagator_.PropagateWithStm(odefunc_, t0, tf, x0, dt, J);
    stm = J;
    return xf;
  }

  Vec6 NumericalOrbitDynamics::PropagateR(Vec6 &x, Real t0, Real tf) {
    Real dt_prop = (dt_ > 0.0) ? dt_ : (tf - t0) / 10;
    return PropagateR(x, t0, tf, dt_prop);
  }

  Vec6 NumericalOrbitDynamics::PropagateWithStmR(Vec6 &x, Real t0, Real tf, Mat6d &stm) {
    Real dt_prop = (dt_ > 0.0) ? dt_ : (tf - t0) / 10;
    return PropagateWithStmR(x, t0, tf, dt_prop, stm);
  }

  // ****************************************************************************
  // CartesianTwoBodyDynamics
  // ****************************************************************************

  CartesianTwoBodyDynamics::CartesianTwoBodyDynamics(double GM, std::string integratorType)
      : mu_(GM),
        NumericalOrbitDynamics(std::bind(&CartesianTwoBodyDynamics::ComputeRates, this,
                                         std::placeholders::_1, std::placeholders::_2),
                               OrbitStateRepres::CARTESIAN, integratorType) {};

  VecX CartesianTwoBodyDynamics::ComputeRates(Real t, const VecX &x) const {
    Vec6 dxdt;
    Vec3 r = x.head(3);
    Vec3 v = x.tail(3);
    Real r_norm = r.norm();

    dxdt.head(3) = v;
    dxdt.tail(3) = -mu_ * r / pow(r_norm, 3);

    return dxdt;
  }

  // ****************************************************************************
  // MoonFixedDynamics
  // ****************************************************************************

  MoonFixedDynamics::MoonFixedDynamics(double GM, std::string integratorType)
      : mu_(GM),
        NumericalOrbitDynamics(std::bind(&MoonFixedDynamics::ComputeRates, this,
                                         std::placeholders::_1, std::placeholders::_2),
                               OrbitStateRepres::CARTESIAN, integratorType) {};

  VecX MoonFixedDynamics::ComputeRates(Real t, const VecX &x) const {
    Vec6 dxdt;
    Vec3 r = x.head(3);
    Vec3 v = x.tail(3);
    Real r_norm = r.norm();

    dxdt.head(3) = Vec3::Zero();
    dxdt.tail(3) = Vec3::Zero();

    return dxdt;
  }

  // ****************************************************************************
  // J2CartTwoBodyDynamics
  // ****************************************************************************

  J2CartTwoBodyDynamics::J2CartTwoBodyDynamics(double GM, double J2, double Rbody,
                                               std::string integratorType)
      : mu_(GM),
        J2_(J2),
        Rbody_(Rbody),
        NumericalOrbitDynamics(std::bind(&J2CartTwoBodyDynamics::ComputeRates, this,
                                         std::placeholders::_1, std::placeholders::_2),
                               OrbitStateRepres::CARTESIAN, integratorType) {};

  VecX J2CartTwoBodyDynamics::ComputeRates(Real t, const VecX &x) const {
    VecX acc(6);

    Vec3 r = x.head(3);
    Vec3 v = x.tail(3);
    Real r_norm = r.norm();

    acc.head(3) = v;
    acc.tail(3) = -mu_ * r / pow(r_norm, 3);

    Vec3 aJ2;
    Real aux1 = -3.0 / 2.0 * mu_ * J2_ * pow(Rbody_, 2.0) / pow(r_norm, 4.0);
    Real aux2 = 5.0 * pow(r(2) / Rbody_, 2.0);
    aJ2(0) = aux1 * (1.0 - aux2) * r(0);
    aJ2(1) = aux1 * (1.0 - aux2) * r(1);
    aJ2(2) = aux1 * (3.0 - aux2) * r(2);

    acc.tail(3) += aJ2;

    return acc;
  }

  // ****************************************************************************
  // J2KeplerianDynamics
  // ****************************************************************************

  J2KeplerianDynamics::J2KeplerianDynamics(double mu_in, double J2_in, double Rbody_in,
                                           std::string integratorType)
      : GM(mu_in),
        J2(J2_in),
        Rbody(Rbody_in),
        NumericalOrbitDynamics(std::bind(&J2KeplerianDynamics::ComputeRates, this,
                                         std::placeholders::_1, std::placeholders::_2),
                               OrbitStateRepres::CARTESIAN, integratorType) {};

  VecX J2KeplerianDynamics::ComputeRates(Real t, const VecX &x) const {
    Real p = x(0) * (1.0 - x(1) * x(1));
    Real n = sqrt(GM / pow(x(0), 3.0));
    Real eta = sqrt(1.0 - x(1) * x(1));

    VecX coeDot(6);
    coeDot(0) = 0.0;
    coeDot(1) = 0.0;
    coeDot(2) = 0.0;
    coeDot(3) = -3.0 / 2.0 * J2 * pow(Rbody / p, 2.0) * n * cos(x(2));
    coeDot(4) = 3.0 / 4.0 * J2 * pow(Rbody / p, 2.0) * n * (5.0 * pow(cos(x(2)), 2.0) - 1.0);
    coeDot(5)
        = n + 3.0 / 4.0 * J2 * pow(Rbody / p, 2.0) * n * eta * (3.0 * pow(cos(x(2)), 2.0) - 1.0);
    return coeDot;
  }

  // ****************************************************************************
  // MoonMeanDynamics
  // ****************************************************************************

  MoonMeanDynamics::MoonMeanDynamics(std::string integratorType)
      : NumericalOrbitDynamics(std::bind(&MoonMeanDynamics::ComputeRates, this,
                                         std::placeholders::_1, std::placeholders::_2),
                               OrbitStateRepres::CARTESIAN, integratorType) {};

  VecX MoonMeanDynamics::ComputeRates(Real t, const VecX &x) const {
    Real a = x(0);
    Real e = x(1);
    Real i = x(2);
    Real O = x(3);
    Real w = x(4);
    Real M = x(5);

    VecX coeDot(6);
    coeDot(0) = 0;
    coeDot(1) = (15 * k * pow(n3, 2) * pow(a, 3.0 / 2.0)) / (8 * sqrt(GM_MOON)) * e
                * sqrt(1 - pow(e, 2)) * pow(sin(i), 2) * sin(2 * w);
    coeDot(2) = -(15 * k * pow(n3, 2) * pow(a, 3.0 / 2.0)) / (16 * sqrt(GM_MOON)) * pow(e, 2)
                / sqrt(1 - pow(e, 2)) * sin(2 * i) * sin(2 * w);
    coeDot(3) = -(3 * J2 * sqrt(GM_MOON) * pow(R_MOON, 2))
                    / (2 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 2)) * cos(i)
                + (3 * k * pow(n3, 2) * pow(a, 3.0 / 2.0))
                      / (8 * sqrt(GM_MOON) * sqrt(1 - pow(e, 2)))
                      * (5 * pow(e, 2) * cos(2 * w) - 3 * pow(e, 2) - 2) * cos(i);
    coeDot(4) = (3 * J2 * sqrt(GM_MOON) * pow(R_MOON, 2))
                    / (4 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 2)) * (5 * pow(cos(i), 2) - 1)
                + (3 * k * pow(n3, 2) * pow(a, 3.0 / 2.0))
                      / (8 * sqrt(GM_MOON) * sqrt(1 - pow(e, 2)))
                      * ((5 * pow(cos(i), 2) - 1 + pow(e, 2))
                         + 5 * (1 - pow(e, 2) - pow(cos(i), 2)) * cos(2 * w));
    coeDot(5) = sqrt(GM_MOON) / pow(a, 3.0 / 2.0)
                + (3 * J2 * sqrt(GM_MOON) * pow(R_MOON, 2))
                      / (4 * pow(a, 7.0 / 2.0) * pow(1 - pow(e, 2), 3.0 / 2.0))
                      * (3 * pow(cos(i), 2) - 1)
                - (k * pow(n3, 2) * pow(a, 3.0 / 2.0)) / (8 * sqrt(GM_MOON))
                      * ((3 * pow(e, 2) + 7) * (3 * pow(cos(i), 2) - 1)
                         + 15 * (1 + pow(e, 2)) * pow(sin(i), 2) * cos(2 * w));
    return coeDot;
  }
};  // namespace lupnt
