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

using ODE = std::function<VecX(const real, const VecX &)>;

class IIntegrator {
 public:
  virtual VecX Step(const ODE f, const real t, const VecX x, const real dt) = 0;
  virtual ~IIntegrator() {};
};

class RK4 : public IIntegrator {
 public:
  VecX Step(const ODE f, const real t, const VecX x, const real dt);
};

class RK8 : public IIntegrator {
 public:
  VecX Step(const ODE f, const real t, const VecX x, const real dt);
};

}  // namespace lupnt