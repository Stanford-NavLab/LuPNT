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
  // Dynamics
  auto dyn_earth_tb = std::make_shared<CartesianTwoBodyDynamics>(MU_EARTH);
  auto dyn_moon_tb = std::make_shared<CartesianTwoBodyDynamics>(MU_MOON);
  auto dyn_moon_nb = std::make_shared<NBodyDynamics>();

  auto earth = Body::Earth();
  auto moon = Body::Moon(5, 5);
  dyn_moon_nb->AddBody(moon);
  dyn_moon_nb->AddBody(earth);
  dyn_moon_nb->SetCentralBody(moon);

  auto dyn_clk = ClockDynamics(ClockModel::kMicrosemiCsac);

  // GPS constellation
  auto channel = std::make_shared<GnssChannel>();
  auto gps_const = GnssConstellation();
  gps_const.SetChannel(channel);
  gps_const.SetDynamics(dyn_earth_tb);
  gps_const.LoadTleFile("gps");

  // Time
  // double epoch0 = sp::StringToTAI("2001/04/06 07:51:28.788 UTC").val();
  double epoch0 = gps_const.GetEpoch();
  std::string epoch_string = sp::TAItoStringUTC(epoch0, 3);
  std::cout << "Initial Epoch: " << epoch_string << std::endl;

  double t0 = 0;
  double tf = t0 + 1.0 * SECS_PER_HOUR;  // 6 hours
  double dt = 1.0;                       // Integration time step [s]
  double Dt = 10.0;                      // Propagation time step [s]
  double print_every = 0.1 * SECS_PER_HOUR;
  double save_every = Dt;

  // Moon spacecraft
  real a = 6541.4;
  real e = 0.6;
  real i = 65.5 * RAD_PER_DEG;
  real Omega = 0.0 * RAD_PER_DEG;
  real w = 90.0 * RAD_PER_DEG;
  real M = 0.0 * RAD_PER_DEG;
  ClassicalOE coe_moon({a, e, i, Omega, w, M}, CoordSystem::MI);
  auto cart_state_moon = std::make_shared<CartesianOrbitState>(
      ClassicalToCartesian(coe_moon, MU_MOON));

  Vector2 clock_vec{0.0, 0.1};  // [s, s/s]
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

  // Estimation
  int state_size = 8;              // Pos, vel, bias, drift [km, km/s, s, s/s]
  double pos_err = 1.0;            // Position error [km]
  double vel_err = 1e-3;           // Velocity error [km/s]
  double clk_bias_err = 1e-6;      // Clock bias error [s]
  double clk_drift_err = 1e-9;     // Clock drift error [s/s]
  double sigma_acc = 2e-6;         // Process noise [km/s^2]
  double sigma_range = 5e-3;       // Range measurement noise [km]
  double sigma_range_rate = 1e-6;  // Range rate measurement noise [km/s]

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

  FilterMeasurementFunction meas_func_pos_clk =
      [moon_sat, receiver, state_size, sigma_range](
          const VectorX& x, MatrixXd& H, MatrixXd& R) -> VectorX {
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

    VectorX z = meas.GetPseudorange(epoch, x.head(6), x.tail(2), H);
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

  // Output
  auto data_history = std::make_shared<DataHistory>();
  auto output_path = std::filesystem::current_path() / "output" / "ExampleEKF";
  FileWriter writer(output_path, true);

  // Main loop
  double t = t0;
  double epoch = epoch0;
  std::cout << "Run Simulation" << std::endl;
  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Time  | Pos Err | Vel Err | Clk Bias Err" << std::endl;
  std::cout << "-----------------------------------------" << std::endl;

  while (t < tf) {
    t += Dt;
    epoch += Dt;

    // Propagate True State
    moon_sat->Propagate(epoch);
    gps_const.Propagate(epoch);

    // Get True Measurement
    std::string freq = "L1";
    double epoch = moon_sat->GetEpoch().val();
    auto measall = receiver->GetMeasurement(epoch);
    auto meas = measall.ExtractSignal("L1");
    int meas_size = meas.GetNumMeasurements();
    auto CN0 = meas.GetCN0();
    auto z_true = meas.GetPseudorange();

    // Update EKF
    ekf.Step(t + dt, z_true);

    // Add Data
    AddStateEstimationData(data_history, moon_sat, &ekf, &gps_const, &meas, t,
                           epoch);

    // Print progress
    if (fmod(t, print_every) < 1e-3) {
      PrintProgress(t, moon_sat, &ekf);
    }
  }

  // Plots -----------------------------------------------
  // Plot3DTrajectory(data_history, "true");
  // Plot3DTrajectory(data_history, "est");

  // Plot estimate for all states
  PlotState(data_history, "true");
  PlotState(data_history, "est");

  // Write data ------------------------------------------
  std::cout << "Simulation finished, saving data..." << std::endl;
  writer.WriteData(*data_history);
  std::cout << "Data saved to " << output_path << std::endl;
  return 0;
}
