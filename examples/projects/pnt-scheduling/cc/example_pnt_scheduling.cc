#include <lupnt/lupnt.h>

#include <chrono>
#include <iostream>

using namespace lupnt;
using namespace std::chrono;

int temp() {
  // Epoch
  auto epoch_0 = SpiceInterface::StringToTAI("2025/10/02 00:00:00.000 UTC");

  // Orbital elements
  real sma = 5740;             // [km]  a, Semi-major axis
  real ecc = 0.58;             // [-]   e, Eccentricity
  real inc = deg2rad(54.856);  // [rad] i, Inclination
  real raan = deg2rad(0);      // [rad] W, Right ascension of the ascending node
  real aop = deg2rad(86.322);  // [rad] w, Argument of periapsis
  real ma = deg2rad(0);        // [rad] M, Mean anomaly

  Vector6 coe_sat_OP{sma, ecc, inc, raan, aop, ma};
  Vector6 rv0_sat_OP = ClassicalToCartesian(coe_sat_OP, MU_MOON);
  Vector6 rv0_sat_mi = CoordConverter::Convert(
      epoch_0, rv0_sat_OP, CoordSystem::OP, CoordSystem::MI);

  // Time
  real T = 2 * PI * sqrt(pow(sma, 3.) / MU_MOON);
  real Dt = 20 * SECS_PER_MINUTE;
  real dt = 5 * SECS_PER_MINUTE;
  real tf = 5 * SECS_PER_DAY;
  int N_t = int(tf / Dt);
  VectorX tspan = VectorX::LinSpaced(N_t, 0, tf);
  VectorX epochs = epoch_0 + tspan.array();

  // Dynamics (real)
  NBodyDynamics dynamics;
  dynamics.SetPrimaryBody(Body::Moon());
  dynamics.AddBody(Body::Earth());
  dynamics.SetTimeStep(dt);

  // Propagate
  auto start = high_resolution_clock::now();
  auto rv_moon_sat_mi = dynamics.Propagate(rv0_sat_mi, epoch_0, epochs);
  auto end = high_resolution_clock::now();
  auto duration = duration_cast<seconds>(end - start);
  std::cout << "Time elapsed (real): " << duration.count() << " seconds"
            << std::endl;
}

int main() {
  Vector3 x{1e5, 1e4, 1e3};
  Vector3 y{10, 20, 30};
  std::cout << decimal2dB(x) << std::endl;
  std::cout << dB2decimal(y) << std::endl;
  return 0;
}