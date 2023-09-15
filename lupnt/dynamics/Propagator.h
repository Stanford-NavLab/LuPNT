/**
 * @file Propagator.h
 * @author Stanford NAV LAB
 * @brief Numerical Propagator
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <autodiff/forward/real.hpp>
#include <iostream>
#include <memory>
#include <string>

#include "lupnt/numerics/Integrator.h"

namespace ad = autodiff;

namespace LPT {
class NumericalPropagator {
 public:
  std::unique_ptr<IIntegrator> integrator;  // integrator type

  NumericalPropagator();
  NumericalPropagator(std::string integratorType);

  ad::VectorXreal Propagate(ODE odefunc, ad::real t0, ad::real tf,
                            ad::VectorXreal x0, ad::real dt);
  ad::VectorXreal PropagateWithStm(ODE odefunc, ad::real t0, ad::real tf,
                                   ad::VectorXreal x0, ad::real dt,
                                   Eigen::MatrixXd &J);
};

}  // namespace LPT
