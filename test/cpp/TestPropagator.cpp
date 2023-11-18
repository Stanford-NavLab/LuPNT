#include <gtest/gtest.h>
#include <lupnt/dynamics/propagator.h>

#include "TestUtils.h"
using namespace lupnt;

namespace {

double MU_MOON = 4902.800066;  // [km^3/s^2]

ODE TwoBodyODE = [](const ad::real t, const ad::VectorXreal& x) {
  ad::VectorXreal acc(6);

  ad::Vector3real r = x.head(3);
  ad::Vector3real v = x.tail(3);
  ad::real r_norm = r.norm();

  acc.head(3) = v;
  acc.tail(3) = -MU_MOON * r / pow(r_norm, 3);

  return acc;
};

TEST(NumericalPropagator, PropagateWithStmTwoBodyTest) {
  NumericalPropagator propagator("RK4");

  ad::real t0 = 0;
  ad::VectorXreal x0(6);
  x0 << 1.60218e-13, 2616.56, 0.0, -0.718032, 4.39668e-17, 1.57558;

  ad::VectorXreal x0_num = x0;

  ad::real dt = 1.0;
  ad::real Dt = 10.0;
  Eigen::MatrixXd J;

  for (int i = 0; i < 100; i++) {
    ad::VectorXreal xEnd =
        propagator.PropagateWithStm(TwoBodyODE, t0, t0 + Dt, x0, dt, J);
  }
}

}  // namespace
