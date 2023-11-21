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

namespace ad = autodiff;

namespace lupnt {

using ODE = std::function<VectorXreal(const real, const VectorXreal &)>;

class IIntegrator {
 public:
  virtual VectorXreal Step(const ODE &f, const real t, const VectorXreal &x,
                           const real dt) = 0;
  virtual ~IIntegrator(){};
};

class RK4 : public IIntegrator {
 public:
  VectorXreal Step(const ODE &f, const real t, const VectorXreal &x,
                   const real dt);
};

class RK8 : public IIntegrator {
 public:
  VectorXreal Step(const ODE &f, const real t, const VectorXreal &x,
                   const real dt);
};

}  // namespace lupnt