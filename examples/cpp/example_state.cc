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
  clock.SetValue(0.0, 0);
  clock.SetValue(0.1, 1);

  // Set the state
  real a = 6541.4;
  real e = 0.6;
  real i = 65.5 * RAD_PER_DEG;
  real Omega = 0.0 * RAD_PER_DEG;
  real w = 90.0 * RAD_PER_DEG;
  real M = 0.0 * RAD_PER_DEG;
  ClassicalOE coe({a, e, i, Omega, w, M});
  auto cart = ClassicalToCartesian(coe, MU_MOON);

  // dynamics model
  auto dyn_moon_tb = CartesianTwoBodyDynamics(MU_MOON);
  dyn_moon_tb.SetDt(0.5);
  auto dyn_clk = ClockDynamics(ClockModel::kMicrosemiCsac);

  // Joint state
  JointState joint_state;
  joint_state.PushBackStateAndDynamics(&cart, &dyn_moon_tb);
  joint_state.PushBackStateAndDynamics(&clock, &dyn_clk);

  VectorX state_vec = joint_state.GetJointStateValue();
  std::cout << "State: " << std::endl << state_vec << std::endl;

  // Test propagation
  DynamicsFunction joint_dynamics = joint_state.GetDynamicsFunction();
  MatrixXd Phi;
  real t_start = 0.0;
  real t_end = 60.0;
  VectorX prop_state = joint_dynamics(state_vec, t_start, t_end, Phi);

  std::cout << "Propagated State: " << std::endl << prop_state << std::endl;
  std::cout << "State Transition Matrix: " << std::endl << Phi << std::endl;
}