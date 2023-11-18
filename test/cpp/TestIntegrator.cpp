#include <gtest/gtest.h>
#include <lupnt/numerics/integrator.h>

using namespace lupnt;

namespace {

double OMEGA = 2.0 * M_PI;  // angular frequency
double EPS = 1e-5;          // tolerance

ODE HarmonicOscillator = [](const ad::real t, const ad::VectorXreal& x) {
  ad::VectorXreal dxdt(2);
  dxdt[0] = x[1];
  dxdt[1] = -OMEGA * OMEGA * x[0];
  return dxdt;
};

TEST(Integrator, RK4) {
  RK4 integrator;

  ad::real t = 0;
  ad::VectorXreal x(2);
  x[0] = 1;            // initial position
  x[1] = 0;            // initial velocity
  ad::real dt = 0.01;  // time step

  // Simulate for one period
  for (int i = 0; i < (1.0 / dt); i++) {
    x = integrator.Step(HarmonicOscillator, t, x, dt);
    t += dt;
  }

  // After one period, the oscillator should be back to the initial position.
  EXPECT_NEAR(x(0).val(), 1.0, EPS);
  EXPECT_NEAR(x(1).val(), 0.0, EPS);
}

TEST(Integrator, RK8) {
  RK8 integrator;

  ad::real t = 0;
  ad::VectorXreal x(2);
  x[0] = 1;            // initial position
  x[1] = 0;            // initial velocity
  ad::real dt = 0.01;  // time step

  // Simulate for one period
  for (int i = 0; i < (1.0 / dt); i++) {
    x = integrator.Step(HarmonicOscillator, t, x, dt);
    t += dt;
  }

  // After one period, the oscillator should be back to the initial position.
  EXPECT_NEAR(x(0).val(), 1.0, EPS);
  EXPECT_NEAR(x(1).val(), 0.0, EPS);
}

}  // namespace
