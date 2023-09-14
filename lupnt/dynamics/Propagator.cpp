#include "Propagator.h"

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

namespace ad = autodiff;

namespace LPT {

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

ad::VectorXreal NumericalPropagator::Propagate(ODE odefunc, ad::real t0,
                                               ad::real tf, ad::VectorXreal x0,
                                               ad::real dt) {
  ad::VectorXreal x = x0;
  ad::real t = t0;
  ad::real step;
  while (t <= tf) {
    step = std::min(dt, tf - t);
    x = integrator->Step(odefunc, t, x, step);
    t += dt;
  }
  return x;
};

ad::VectorXreal NumericalPropagator::PropagateWithStm(ODE odefunc, ad::real t0,
                                                      ad::real tf,
                                                      ad::VectorXreal x0,
                                                      ad::real dt,
                                                      Eigen::MatrixXd &J) {
  auto func = [=, this](ad::VectorXreal &x) {
    return Propagate(odefunc, t0, tf, x, dt);
  };

  // decouple x0 from previous relations by reinitializing it
  for (int i = 0; i < x0.size(); i++) {
    x0(i) = double(x0(i));
  }
  ad::VectorXreal xf;
  J = ad::jacobian(func, wrt(x0), at(x0), xf);
  return xf;
};

};  // namespace LPT