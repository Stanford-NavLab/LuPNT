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

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

namespace ad = autodiff;

namespace lupnt {

typedef std::function<ad::VectorXreal(const ad::real, const ad::VectorXreal &)>
    ODE;

class IIntegrator {
 public:
  virtual ad::VectorXreal Step(const ODE &f, const ad::real t,
                               const ad::VectorXreal &x, const ad::real dt) = 0;
  virtual ~IIntegrator(){};
};

class RK4 : public IIntegrator {
 public:
  ad::VectorXreal Step(const ODE &f, const ad::real t, const ad::VectorXreal &x,
                       const ad::real dt);
};

class RK8 : public IIntegrator {
 public:
  ad::VectorXreal Step(const ODE &f, const ad::real t, const ad::VectorXreal &x,
                       const ad::real dt);
};

}  // namespace lupnt