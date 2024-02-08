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

#include <cmath>

using namespace lupnt;

// Define Measurement and Process Noise Functions
int state_size = 8;
int meas_size = 4;

MatrixXd SetInitialCovariance(double pos_sigma, double vel_sigma,
                              double clock_bias_sigma,
                              double clock_drift_sigma) {
  MatrixXd P0 = MatrixXd::Identity(state_size, state_size);
  P0.block(0, 0, 3, 3) = Matrix3d::Identity() * pos_sigma;
  P0.block(3, 3, 3, 3) = Matrix3d::Identity() * vel_sigma;
  P0(6, 6) = clock_bias_sigma;
  P0(7, 7) = clock_drift_sigma;
  return P0;
}

FilterProcessNoiseFunction proc_noise_func = [](const VectorX& x, real t_curr,
                                                real t_end) {
  int clock_index = 6;
  double dt = (t_end - t_curr).val();
  double sigma_acc = 1e-13;

  MatrixXd Q = MatrixXd::Zero(state_size, state_size);

  Matrix6d Q_rv = Matrix6d::Zero();
  for (int i = 0; i < 3; i++) {
    Q_rv(i, i) = pow(dt, 3) / 3.0 * pow(sigma_acc, 2);
    Q_rv(i + 3, i + 3) = dt * pow(sigma_acc, 2);
    Q_rv(i, i + 3) = pow(dt, 2) / 2.0 * pow(sigma_acc, 2);
    Q_rv(i + 3, i) = pow(dt, 2) / 2.0 * pow(sigma_acc, 2);
  }
  Matrix2d Q_clk = GetClockProcessNoise(ClockModel::kMicrosemiCsac, dt);

  Q.block(0, 0, 6, 6) = Q_rv;
  Q.block(6, 6, 2, 2) = Q_clk;

  return Q;
};

FilterMeasurementFunction meas_func_pos_clk = [](const VectorX& x, MatrixXd& H,
                                                 MatrixXd& R) {
  int clock_index = 6;

  VectorX y = VectorX::Zero(meas_size);
  H.resize(meas_size, state_size);
  R.resize(meas_size, meas_size);

  Vector3 pos = x.segment(0, 3);
  Vector3 vel = x.segment(3, 3);
  H = MatrixXd::Zero(meas_size, state_size);
  H.block(0, 0, 3, 3) = Matrix3d::Identity();
  H(3, clock_index) = 1.0;
  R = Matrix4d::Identity() * 1e-3;

  y.segment(0, 3) = pos;
  y(3) = x(clock_index);

  return y;
};

// Main
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
  dyn_moon_tb.SetTimeStep(0.5);
  auto dyn_clk = ClockDynamics(ClockModel::kMicrosemiCsac);

  // Joint state
  JointState joint_state;
  joint_state.PushBackStateAndDynamics(&cart, &dyn_moon_tb);
  joint_state.PushBackStateAndDynamics(&clock, &dyn_clk);

  VectorX state_vec = joint_state.GetJointStateValue();
  std::cout << "Initial State: " << std::endl << state_vec << std::endl;

  // Test propagation
  FilterDynamicsFunction joint_dynamics =
      joint_state.GetFilterDynamicsFunction();
  MatrixXd Phi;
  real t_start = 0.0;
  real t_end = 60.0;
  VectorX prop_state = joint_dynamics(state_vec, t_start, t_end, Phi);

  std::cout << " " << std::endl;
  std::cout << "Propagated State (60 sec): " << std::endl
            << prop_state << std::endl;
  std::cout << " " << std::endl;
  std::cout << "State Transition Matrix: " << std::endl << Phi << std::endl;

  // Initialize EKF
  EKF ekf;
  ekf.SetDynamicsFunction(joint_dynamics);
  ekf.SetMeasurementFunction(meas_func_pos_clk);
  ekf.SetProcessNoiseFunction(proc_noise_func);

  // Simulation Time
  t_start = 0.0;
  t_end = 100.0;
  real dt = 2.0;

  // Dummy Variables for Storage
  MatrixXd Phat(state_size, state_size);
  MatrixXd Phi_t(state_size, state_size);
  MatrixXd H(meas_size, state_size);
  MatrixXd R(meas_size, meas_size);

  // Initialize
  MatrixXd P0 = SetInitialCovariance(1e-3, 1e-6, 1e-9, 1e-12);  // Initial error
  VectorX x_true(state_size), x_est(state_size);
  x_true = state_vec;
  x_est = x_true + SampleMVN(Vector8d::Zero(), P0, 1).col(0);
  ekf.Initialize(x_est, P0);

  MatrixXd Q = ekf.process_noise_(x_true, 0, dt);
  std::cout << " " << std::endl;
  std::cout << "Process Noise: " << std::endl << Q << std::endl;

  // Run EKF Simulation
  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Time  | Pos Err | Vel Err | Clk Bias Err" << std::endl;
  std::cout << "-----------------------------------------" << std::endl;

  for (real t = t_start; t < t_end; t += dt) {
    // Propagate True State
    x_true = joint_dynamics(x_true, t, t + dt, Phi_t);
    x_true += SampleMVN(Vector8d::Zero(), Q, 1).col(0);

    // Get True Measurement
    VectorX meas = ekf.measurement_(x_true, H, R);
    meas += SampleMVN(Vector4d::Zero(), R, 1).col(0);

    // Update EKF
    ekf.Step(t + dt, meas);

    // Print Result every 10 seconds
    if (std::fmod(t.val(), 10.0) != 0.0) {
      continue;
    } else {
      x_est = ekf.GetUpdatedStateEstimate(Phat);
      double x_pos_err =
          (x_true.segment(0, 3) - x_est.segment(0, 3)).norm().val();
      double x_vel_err =
          (x_true.segment(3, 3) - x_est.segment(3, 3)).norm().val();
      double x_clk_bias_err = abs((x_true(6) - x_est(6)).val());

      std::cout.precision(3);
      std::cout << t << " | " << x_pos_err << " | " << x_vel_err << " | "
                << x_clk_bias_err << std::endl;
    }
  }

  // Print final covariance
  std::cout << " " << std::endl;
  std::cout << "Final Covariance: " << std::endl << Phat << std::endl;
}