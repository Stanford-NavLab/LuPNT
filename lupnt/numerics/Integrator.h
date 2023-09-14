#pragma once

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

namespace ad = autodiff;

namespace LPT {

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

}  // namespace LPT