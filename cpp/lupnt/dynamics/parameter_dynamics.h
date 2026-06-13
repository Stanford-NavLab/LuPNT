/**
 * @file parameter_dynamics.h
 * @author Stanford NAV LAB
 * @brief  Parameter dynamics classes
 * @version 0.1
 * @date 2025-11-07
 */

#pragma once

#include "lupnt/dynamics/dynamics.h"
#include "lupnt/states/params.h"

namespace lupnt {

  /**
   * @class ParameterDynamics
   * @brief Class for parameter dynamics.
   *
   * This class provides methods for propagating parameter states, which are assumed to be constant.
   */
  class ParameterDynamics : public Dynamics {
  public:
    ParameterDynamics() = default;
    ParameterDynamics(Config& config) : Dynamics(config) {}

    State Propagate(const State& x0, Real t0, Real tf, const State* u = nullptr) override {
      (void)t0;
      (void)tf;
      (void)u;
      return x0;  // Parameters are constant
    }

    State Propagate(const State& x0, Real t0, Real tf, const State* u, MatXd* stm) override {
      (void)t0;
      (void)tf;
      (void)u;
      if (stm != nullptr) {
        *stm = MatXd::Identity(x0.size(), x0.size());
      }
      return x0;  // Parameters are constant
    }

    StateType GetStateType() const override { return "Params"; }
  };

}  // namespace lupnt
