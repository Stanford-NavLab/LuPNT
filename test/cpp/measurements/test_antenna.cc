#include <lupnt/lupnt.h>

#include "../utils.cc"

double ABS_TOL = 1e-6;
double REL_TOL = 1e-3;

TEST_CASE("Antenna") {
  using namespace lupnt;
  Antenna ant("Block-IIR-M_ACE");
  RequireNear(ant.ComputeGain(74 * RAD, 358 * RAD), -18.68, ABS_TOL);
  RequireNear(ant.ComputeGain(74 * RAD, 359 * RAD), -18.73, ABS_TOL);
  RequireNear(ant.ComputeGain(74 * RAD, 358.5 * RAD), 0.5 * (-18.68 - 18.73), ABS_TOL);
  RequireNear(ant.ComputeGain(74 * RAD, 359.5 * RAD), 0.5 * (-20.45 - 18.73), ABS_TOL);
  RequireNear(ant.ComputeGain(74.5 * RAD, 359 * RAD), 0.5 * (-17.48 - 18.73), ABS_TOL);
}
