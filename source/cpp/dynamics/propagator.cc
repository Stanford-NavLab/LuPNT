/**
 * @file Propagator.cpp
 * @author Stanford NAV LAB
 * @brief Propagator class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "lupnt/dynamics/propagator.h"

#include <algorithm>

namespace lupnt {

  NumericalPropagator::NumericalPropagator(IntegratorType integ, IntegratorParams params) {
    if (integ == IntegratorType::RK4)
      integrator = MakePtr<RK4>();
    else if (integ == IntegratorType::RK8)
      integrator = MakePtr<RK8>();
    else if (integ == IntegratorType::RKF45)
      integrator = MakePtr<RKF45>(params);
    else
      throw std::invalid_argument("Invalid Integrator Type");
  };

  VecX NumericalPropagator::Propagate(const ODE &odefunc, Real t0, Real tf, const VecX &x0,
                                      Real dt) {
    assert(dt > 0 && "Invalid time step");

    VecX x = x0;
    Real t = t0;
    while (t < tf) {
      dt = std::min(dt, tf - t);
      // store the previous step (deep copy)
      Real prev_dt = dt;
      x = integrator->Step(odefunc, t, x, dt);  // update x and step
      t += prev_dt;
      // std::cout << "t: " << t << std::endl;
      if (log_history_) {
        t_history_.push_back(t);
        x_history_.push_back(x);
      }
    }
    return x;
  };

  VecX NumericalPropagator::Propagate(const ODE &odefunc, Real t0, Real tf, const VecX &x0, Real dt,
                                      MatXd *J) {
    auto func = [=, this](const VecX &x) { return Propagate(odefunc, t0, tf, x, dt); };

    VecX xf;
    if (J != nullptr) {
      // decouple x0 from previous relations by reinitializing it
      VecX x0_tmp = x0.cast<double>();
      *J = jacobian(func, wrt(x0_tmp), at(x0_tmp), xf);
    } else {
      xf = func(x0);
    }
    return xf;
  };

};  // namespace lupnt
