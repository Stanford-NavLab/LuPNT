/**
 * @file ExampleEKF.cpp
 * @author Stanford NAV Lab
 * @brief Example of Extended Kalman Filter for orbit determination
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

// lupnt includes
#include "example_utils.cc"

using namespace lupnt;
namespace sp = SpiceInterface;

// Util Functions

int main() {
  /**********************************************
   * Simulation Parameters
   *********************************************/
  // time
  double t0 = 0;
  double tf = t0 + 1.0 * SECS_PER_HOUR;  // 6 hours
  double dt = 1.0;                       // Integration time step [s]
  double Dt = 10.0;                      // Propagation time step [s]
  double print_every = 0.1 * SECS_PER_HOUR;
  double save_every = Dt;

  // for debug
  tf = t0 + 3.0 * SECS_PER_HOUR;
  print_every = 600;

  // Initial State
  real a = 6541.4;
  real e = 0.6;
  real i = 65.5 * RAD_PER_DEG;
  real Omega = 0.0 * RAD_PER_DEG;
  real w = 90.0 * RAD_PER_DEG;
  real M = 0.0 * RAD_PER_DEG;
  real clk_bias = 0.0;
  real clk_drift = 0.1;

  // Estimation
  int state_size = 8;              // Pos, vel, bias, drift [km, km/s, s, s/s]
  double pos_err = 1.0;            // Position error [km]
  double vel_err = 1e-3;           // Velocity error [km/s]
  double clk_bias_err = 1e-6;      // Clock bias error [s]
  double clk_drift_err = 1e-9;     // Clock drift error [s/s]
  double sigma_acc = 2e-6;         // Process noise [km/s^2]
  double sigma_range = 5e-3;       // Range measurement noise [km]
  double sigma_range_rate = 1e-6;  // Range rate measurement noise [km/s]

  // Debug mode
  bool print_debug = false;
  bool no_meas = false;  // set to true to turn off measurements
  bool use_nbody = false;

  /**********************************************
   * Setup
   * ********************************************/

  // Dynamics
  auto dyn_earth_tb = std::make_shared<CartesianTwoBodyDynamics>(MU_EARTH);
  auto dyn_moon_tb = std::make_shared<CartesianTwoBodyDynamics>(MU_MOON);
  auto dyn_moon_nb = std::make_shared<NBodyDynamics>();

  if (use_nbody) {
    auto earth = Body::Earth();
    auto moon = Body::Moon(5, 5);
    dyn_moon_nb->AddBody(moon);
    dyn_moon_nb->AddBody(earth);
    dyn_moon_nb->SetCentralBody(moon);
  } else {
    dyn_moon_nb->SetCentralBody(Body::Moon());
    dyn_moon_nb->AddBody(Body::Moon());
  }

  auto dyn_clk = ClockDynamics(ClockModel::kMicrosemiCsac);

  // GPS constellation
  auto channel = std::make_shared<GnssChannel>();
  auto gps_const = GnssConstellation();
  gps_const.SetChannel(channel);
  gps_const.SetDynamics(dyn_earth_tb);
  gps_const.LoadTleFile("gps");  // example gps file

  // Time
  double epoch0 = gps_const.GetEpoch();
  std::string epoch_string = sp::TAItoStringUTC(epoch0, 3);
  std::cout << "Initial Epoch: " << epoch_string << std::endl;

  // Set dynamics integration time
  dyn_earth_tb->SetDt(dt);
  dyn_moon_tb->SetDt(dt);
  dyn_moon_nb->SetDt(dt);

  // Moon spacecraft
  ClassicalOE coe_moon({a, e, i, Omega, w, M}, CoordSystem::MI);
  auto cart_state_moon = std::make_shared<CartesianOrbitState>(
      ClassicalToCartesian(coe_moon, MU_MOON));

  Vector2 clock_vec{clk_bias, clk_drift};  // [s, s/s]
  ClockState clock_state(clock_vec);

  auto moon_sat = std::make_shared<Spacecraft>();
  auto receiver = std::make_shared<GnssReceiver>("moongpsr");

  moon_sat->AddDevice(receiver);
  moon_sat->SetDynamics(dyn_moon_nb);
  moon_sat->SetClock(clock_state);
  moon_sat->SetOrbitState(cart_state_moon);
  moon_sat->SetEpoch(epoch0);
  moon_sat->SetBodyId(BodyId::MOON);
  moon_sat->SetClockDynamics(dyn_clk);

  receiver->SetAgent(moon_sat);
  receiver->SetChannel(channel);
  channel->AddReceiver(receiver);

  // Initial covariance
  MatrixXd P0 =
      ConstructInitCovariance(pos_err, vel_err, clk_bias_err, clk_drift_err);

  // Joint state and dynamics
  JointState joint_state;
  joint_state.PushBackStateAndDynamics(cart_state_moon.get(),
                                       dyn_moon_tb.get());
  joint_state.PushBackStateAndDynamics(&clock_state, &dyn_clk);

  FilterDynamicsFunction joint_dynamics =
      joint_state.GetFilterDynamicsFunction();

  // Define Measurement function
  FilterMeasurementFunction meas_func_pos_clk =
      [moon_sat, receiver, state_size, sigma_range, no_meas](
          const VectorX x, MatrixXd& H, MatrixXd& R) -> VectorX {
    if (no_meas) {
      return VectorXd::Zero(0);
    }

    // Measurements
    std::string freq = "L1";
    double epoch = moon_sat->GetEpoch().val();
    auto measall = receiver->GetMeasurement(epoch);
    auto meas = measall.ExtractSignal("L1");
    int meas_size = meas.GetNumMeasurements();
    auto CN0 = meas.GetCN0();

    // Predict measurements
    H = MatrixXd::Zero(meas_size, state_size);
    R = MatrixXd::Zero(meas_size, meas_size);

    // create new state
    VectorX x_rv(6), x_clk(2);
    x_rv = x.head(6);
    x_clk = x.tail(2);

    VectorX z = meas.GetPseudorange2(epoch, x_rv, x_clk, H);
    R.diagonal().array() = pow(sigma_range, 2);

    return z;
  };

  // EKF
  EKF ekf;
  ekf.SetDynamicsFunction(joint_dynamics);
  ekf.SetMeasurementFunction(meas_func_pos_clk);
  ekf.SetProcessNoiseFunction(proc_noise_func);
  std::cout << "Initialized EKF" << std::endl;

  VectorX x_est = SampleMVN(joint_state.GetJointStateValue(), P0, 1);
  ekf.Initialize(x_est, P0);

  // Print State
  if (print_debug) {
    std::cout << "Initial true state: "
              << moon_sat->GetStateVector().transpose() << std::endl;
    std::cout << "Initial estimated state: " << x_est.transpose() << std::endl;
  }

  // Output
  auto data_history = std::make_shared<DataHistory>();
  auto output_path = std::filesystem::current_path() / "output" / "ExampleEKF";
  FileWriter writer(output_path, true);

  /*************************
   * Main loop
   *************************/
  real t = t0;
  double epoch = epoch0;
  PrintProgressHeader();

  for (t = t0; t < tf; t += Dt) {
    epoch += Dt;  // first propagate to the next epoch

    // Propagate True State
    moon_sat->Propagate(epoch);
    gps_const.Propagate(epoch);

    // Get True Measurement
    VectorX z_true;
    auto measall = receiver->GetMeasurement(epoch);
    auto meas = measall.ExtractSignal("L1");
    int meas_size = meas.GetNumMeasurements();
    auto CN0 = meas.GetCN0();

    if (!no_meas) {
      z_true = meas.GetPseudorange();
      // Add Noise
      MatrixX R_noise(meas_size, meas_size);
      R_noise = pow(sigma_range, 2) * MatrixX::Identity(meas_size, meas_size);
      z_true += SampleMVN(VectorX::Zero(meas_size), R_noise, 1);
    } else {
      z_true = VectorXd::Zero(0);
    }

    // Update EKF
    ekf.Step(t + Dt, z_true);

    // Add Data
    if (!no_meas) {
      AddStateEstimationData(data_history, moon_sat, &ekf, &gps_const, &meas,
                             t.val(), epoch);
    }

    // Print progress
    if (fmod(t.val(), print_every) < 1e-3) {
      PrintProgress(t.val(), moon_sat, &ekf);
    }

    // print measurement residuals
    if ((print_debug) && (!no_meas)) {
      std::cout << "  Pbar : " << ekf.Pbar_.diagonal().transpose() << std::endl;
      std::cout << "  Q:  " << std::endl << ekf.Q_ << std::endl;
      std::cout << "  Kalman Gain: " << std::endl << ekf.K_ << std::endl;
      std::cout << "  H:  " << std::endl << ekf.H_ << std::endl;
      std::cout << "  R:  " << std::endl << ekf.R_ << std::endl;
      std::cout << "  S:  " << std::endl << ekf.S_ << std::endl;
      std::cout << "  Meas   Residuals: " << ekf.dy_.transpose() << std::endl;
      std::cout << "  Linear Residuals: "
                << (ekf.dy_ - ekf.H_ * ekf.dx_).transpose() << std::endl;
      std::cout << "  dx: " << ekf.dx_.transpose() << std::endl;
      std::cout << "  Phat: " << ekf.P_.diagonal().transpose() << std::endl;
      std::cout << "  " << std::endl;
    }
  }

  // Write data ------------------------------------------
  std::cout << "Simulation finished, saving data..." << std::endl;
  writer.WriteData(*data_history);
  std::cout << "Data saved to " << output_path << std::endl;

  // Plots -----------------------------------------------
  // Plot3DTrajectory(data_history, "true");
  // Plot3DTrajectory(data_history, "est");

  // Plot estimate for all states
  // PlotState(data_history, "true");
  // PlotState(data_history, "est");

  return 0;
}
