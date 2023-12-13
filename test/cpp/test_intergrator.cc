#include <lupnt/numerics/integrator.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace lupnt;
using namespace Catch::Matchers;

namespace {

double OMEGA = 2.0 * M_PI;  // angular frequency
double EPS = 1e-5;          // tolerance

ODE HarmonicOscillator = [](const real t, const VectorX& x) {
  VectorX dxdt(2);
  dxdt[0] = x[1];
  dxdt[1] = -OMEGA * OMEGA * x[0];
  return dxdt;
};

TEST_CASE("Integrator-RK4") {
  RK4 integrator;

  real t = 0;
  VectorX x(2);
  x[0] = 1;        // initial position
  x[1] = 0;        // initial velocity
  real dt = 0.01;  // time step

  // Simulate for one period
  for (int i = 0; i < (1.0 / dt); i++) {
    x = integrator.Step(HarmonicOscillator, t, x, dt);
    t += dt;
  }

  // After one period, the oscillator should be back to the initial position.
  REQUIRE_THAT(x(0).val(), WithinAbs(1.0, EPS));
  REQUIRE_THAT(x(1).val(), WithinAbs(0.0, EPS));
}

TEST_CASE("Integrator-RK8") {
  RK8 integrator;

  real t = 0;
  VectorX x(2);
  x[0] = 1;        // initial position
  x[1] = 0;        // initial velocity
  real dt = 0.01;  // time step

  // Simulate for one period
  for (int i = 0; i < (1.0 / dt); i++) {
    x = integrator.Step(HarmonicOscillator, t, x, dt);
    t += dt;
  }

  // After one period, the oscillator should be back to the initial position.
  REQUIRE_THAT(x(0).val(), WithinAbs(1.0, EPS));
  REQUIRE_THAT(x(1).val(), WithinAbs(0.0, EPS));
}

}  // namespace
