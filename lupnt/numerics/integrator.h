/**
 * @file integrator.h
 * @author Stanford NAV LAB
 * @brief Integrator interfaces
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include "lupnt/core/constants.h"

namespace lupnt {

using ODE = std::function<VecX(const Real, const VecX &)>;

class IIntegrator {
 public:
  virtual VecX Step(const ODE f, const Real t, const VecX x, const Real dt) = 0;
  virtual ~IIntegrator() {};
};

class RK4 : public IIntegrator {
 public:
  VecX Step(const ODE f, const Real t, const VecX x, const Real dt);
};

class RK8 : public IIntegrator {
 public:
  VecX Step(const ODE f, const Real t, const VecX x, const Real dt);
};

}  // namespace lupnt