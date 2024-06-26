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

#include "lupnt/numerics/integrator.h"

namespace lupnt {
class NumericalPropagator {
 public:
  std::unique_ptr<IIntegrator> integrator;  // integrator type

  NumericalPropagator();
  NumericalPropagator(std::string integratorType);

  VecX Propagate(ODE odefunc, real t0, real tf, VecX x0, real dt);
  VecX PropagateWithStm(ODE odefunc, real t0, real tf, VecX x0, real dt,
                        MatXd &J);
};

}  // namespace lupnt
