#include <gtest/gtest.h>
#include <lupnt/dynamics/propagator.h>

#include "TestUtils.h"
using namespace lupnt;

namespace {

double MU_MOON = 4902.800066;  // [km^3/s^2]

ODE TwoBodyODE = [](const real t, const VectorXreal& x) {
  VectorXreal acc(6);

  Vector3real r = x.head(3);
  Vector3real v = x.tail(3);
  real r_norm = r.norm();

  acc.head(3) = v;
  acc.tail(3) = -MU_MOON * r / pow(r_norm, 3);

  return acc;
};

TEST(NumericalPropagator, PropagateWithStmTwoBodyTest) {
  NumericalPropagator propagator("RK4");

  real t0 = 0;
  VectorXreal x0(6);
  x0 << 1.60218e-13, 2616.56, 0.0, -0.718032, 4.39668e-17, 1.57558;

  VectorXreal x0_num = x0;

  real dt = 1.0;
  real Dt = 10.0;
  MatrixXd J;

  for (int i = 0; i < 100; i++) {
    VectorXreal xEnd =
        propagator.PropagateWithStm(TwoBodyODE, t0, t0 + Dt, x0, dt, J);
  }
}

}  // namespace
