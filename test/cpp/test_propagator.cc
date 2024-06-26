#include <lupnt/dynamics/propagator.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "test_utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {

double MU_MOON = 4902.800066;  // [km^3/s^2]

ODE TwoBodyODE = [](const real t, const VecX& x) {
  VecX acc(6);

  Vec3 r = x.head(3);
  Vec3 v = x.tail(3);
  real r_norm = r.norm();

  acc.head(3) = v;
  acc.tail(3) = -MU_MOON * r / pow(r_norm, 3);

  return acc;
};

TEST_CASE("NumericalPropagator", "PropagateWithStmTwoBodyTest") {
  NumericalPropagator propagator("RK4");

  real t0 = 0;
  VecX x0(6);
  x0 << 1.60218e-13, 2616.56, 0.0, -0.718032, 4.39668e-17, 1.57558;

  VecX x0_num = x0;

  real dt = 1.0;
  real Dt = 10.0;
  VecXd J;

  for (int i = 0; i < 100; i++) {
    VecX xEnd = propagator.PropagateWithStm(TwoBodyODE, t0, t0 + Dt, x0, dt, J);
  }
}

}  // namespace
