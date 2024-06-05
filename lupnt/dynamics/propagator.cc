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
#include "propagator.h"

namespace lupnt {

NumericalPropagator::NumericalPropagator() {
  integrator = std::make_unique<RK4>();
};

NumericalPropagator::NumericalPropagator(std::string integratorType) {
  if (integratorType == "RK4")
    integrator = std::make_unique<RK4>();
  else if (integratorType == "RK8")
    integrator = std::make_unique<RK8>();
  else
    throw std::invalid_argument("Invalid Integrator Type");
};

VectorX NumericalPropagator::Propagate(ODE odefunc, real t0, real tf,
                                       VectorX x0, real dt) {
  assert(dt > 0 && "dt must be greater than 0");
  VectorX x = x0;
  real t = t0;
  real step;
  while (t <= tf) {
    step = std::min(dt, tf - t);
    x = integrator->Step(odefunc, t, x, step);
    t += dt;
  }
  return x;
};

VectorX NumericalPropagator::PropagateWithStm(ODE odefunc, real t0, real tf,
                                              VectorX x0, real dt,
                                              MatrixXd &J) {
  auto func = [=, this](VectorX &x) {
    return Propagate(odefunc, t0, tf, x, dt);
  };

  // decouple x0 from previous relations by reinitializing it
  for (int i = 0; i < x0.size(); i++) {
    x0(i) = double(x0(i));
  }
  VectorX xf;
  J = jacobian(func, wrt(x0), at(x0), xf);
  return xf;
};

};  // namespace lupnt