/**
 * @file propagator.h
 * @author Stanford NAV LAB
 * @brief Numerical Propagator
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "lupnt/numerics/integrator.h"

namespace lupnt {
  class NumericalPropagator {
  private:
    std::vector<Real> t_history_;
    std::vector<VecX> x_history_;
    bool log_history_ = false;

  public:
    Ptr<IIntegrator> integrator;  // integrator type

    NumericalPropagator(IntegratorType integ = default_integrator,
                        IntegratorParams params = IntegratorParams());

    VecX Propagate(const ODE &odefunc, Real t0, Real tf, const VecX &x0, Real dt);
    VecX Propagate(const ODE &odefunc, Real t0, Real tf, const VecX &x0, Real dt, MatXd *J);
    void SetLogHistory(bool log_history) { log_history_ = log_history; };
    void GetTimeHistory(std::vector<Real> &t_history) { t_history = t_history_; };
    void GetStateHistory(std::vector<VecX> &x_history) { x_history = x_history_; };
  };

}  // namespace lupnt
