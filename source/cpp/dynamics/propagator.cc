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

  NumericalPropagator::NumericalPropagator() { integrator = std::make_unique<RK4>(); };

  NumericalPropagator::NumericalPropagator(std::string integratorType) {
    IntegratorParams params = IntegratorParams();  // default params
    if (integratorType == "RK4")
      integrator = std::make_unique<RK4>();
    else if (integratorType == "RK8")
      integrator = std::make_unique<RK8>();
    else if (integratorType == "RKF45")
      integrator = std::make_unique<RKF45>(params);
    else
      throw std::invalid_argument("Invalid Integrator Type");
  };

  NumericalPropagator::NumericalPropagator(std::string integratorType, IntegratorParams params) {
    if (integratorType == "RK4")
      integrator = std::make_unique<RK4>();
    else if (integratorType == "RK8")
      integrator = std::make_unique<RK8>();
    else if (integratorType == "RKF45")
      integrator = std::make_unique<RKF45>(params);
    else
      throw std::invalid_argument("Invalid Integrator Type");
  };

  VecX NumericalPropagator::Propagate(ODE odefunc, Real t0, Real tf, VecX x0, Real dt) {
    assert(dt > 0 && "dt must be greater than 0");
    VecX x = x0;
    Real t = t0;
    Real step = dt;
    while (t < tf) {
      step = std::min(step, tf - t);
      // store the previous step (deep copy)
      Real prev_step = step;
      x = integrator->Step(odefunc, t, x, step);  // update x and step
      t += prev_step;
      // std::cout << "t: " << t << std::endl;
      if (log_history_) {
        t_history_.push_back(t);
        x_history_.push_back(x);
      }
    }
    return x;
  };

  VecX NumericalPropagator::PropagateWithStm(ODE odefunc, Real t0, Real tf, VecX x0, Real dt,
                                             MatXd &J) {
    auto func = [=, this](VecX &x) { return Propagate(odefunc, t0, tf, x, dt); };

    // decouple x0 from previous relations by reinitializing it
    for (int i = 0; i < x0.size(); i++) {
      x0(i) = double(x0(i));
    }
    VecX xf;
    J = jacobian(func, wrt(x0), at(x0), xf);
    return xf;
  };

};  // namespace lupnt
