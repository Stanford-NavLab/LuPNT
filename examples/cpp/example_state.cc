/**
 * @file example_state.cc
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2023-12-13
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <lupnt/agents/state_estimation_app.h>
#include <lupnt/core/constants.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/physics/clock.h>
#include <lupnt/physics/orbit_state.h>

using namespace lupnt;

int main() {
  // Create a clock
  ClockState clock(2);
  std::cout << "clock size" << clock.GetSize() << std::endl;

  // Set the state
  real a = 6541.4;
  real e = 0.6;
  real i = 65.5 * RAD_PER_DEG;
  real Omega = 0.0 * RAD_PER_DEG;
  real w = 90.0 * RAD_PER_DEG;
  real M = 0.0 * RAD_PER_DEG;
  ClassicalOE coe({a, e, i, Omega, w, M});

  // dynamics model
  auto dyn_earth_tb = CartesianTwoBodyDynamics(MU_EARTH);
  auto dyn_clk = ClockDynamics(ClockModel::kMicrosemiCsac);

  // Joint state
  JointState joint_state;
  joint_state.PushBackStateAndDynamics(&coe, &dyn_earth_tb);
  joint_state.PushBackStateAndDynamics(&clock, &dyn_clk);

  VectorX state_vec = joint_state.GetJointStateValue();
  std::cout << "State: " << std::endl << state_vec << std::endl;
}