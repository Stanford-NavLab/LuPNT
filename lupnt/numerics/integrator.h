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

using ODE = std::function<VecX(const Real, const VecX&)>;

class IntegratorParams {
 public:
  int max_iter = 20;
  double abstol = 1e-6;
  double reltol = 1e-6;

  IntegratorParams() = default;
  IntegratorParams(int max_iter, double abstol, double reltol)
      : max_iter(max_iter), abstol(abstol), reltol(reltol) {
    CheckIntegratorParams();
  };
  void CheckIntegratorParams();
};

class IIntegrator {
 public:
  virtual VecX Step(const ODE f, const Real t, const VecX x, Real& dt) = 0;
  virtual ~IIntegrator(){};
};

class RK4 : public IIntegrator {
 public:
  VecX Step(const ODE f, const Real t, const VecX x, Real& dt);
};

class RK8 : public IIntegrator {
 public:
  VecX Step(const ODE f, const Real t, const VecX x, Real& dt);
};

class RKF45 : public IIntegrator {
 private:
  IntegratorParams params_;

 public:
  RKF45(IntegratorParams params) : params_(params) {};
  VecX Step(const ODE f, const Real t, const VecX x, Real& dt);
};

}  // namespace lupnt