/**
 * @file example_ekf_elfo.cc
 * @author Stanford NAV Lab
 * @brief   Example of EKF for lunar orbit estimation
 * @version 0.1
 * @date 2024-01-19
 *
 * @copyright Copyright (c) 2024
 *
 */

// lupnt includes
#include <vector>

#include "utils.cc"

using namespace lupnt;
namespace sp = SpiceInterface;

// Util Functions

int main() {
  /**********************************************
   * Simulation Parameters
   *********************************************/
  // Get initial epoch
  auto gps_const = GnssConstellation();
  double epoch0 = gps_const.GetEpoch();

  // Time
  double t0 = epoch0;
  double dt = 1.0;   // Integration time step [s]
  double Dt = 10.0;  // Propagation time step [s]  (= Measurement time step)
  double print_every = 600;
  double save_every = Dt;

  // Simulation seed
  int seed = 0;
  std::srand(seed);

  // Initial State
  real a = 6541.4;
  real e = 0.6;
  real i = 65.5 * RAD_PER_DEG;
  real Omega = 0.0 * RAD_PER_DEG;
  real w = 90.0 * RAD_PER_DEG;
  real M = 0.0 * RAD_PER_DEG;
  real clk_bias = 0.0;
  real clk_drift = 0.1;

  // Set simulation to 1 orbit
  int n_orbit = 1;  // number of orbits to simulate
  real period = 2.0 * M_PI * sqrt(pow(a, 3) / MU_MOON);
  double tf = t0 + n_orbit * period.val();
  int time_step_num = int((tf - t0) / Dt) + 1;
  tf = t0 + (time_step_num - 1) * Dt;

  // Dynamics Model   Todo: Refine this to a more high fidelity model
  int moon_sph_true = 5;  // moon spherical harmonics order in true dynamics
  int moon_sph_est = 0;   // moon spherical harmonics order in filter dynamics
  bool add_earth = true;  // add earth to true and filter dynamics

  // Onboard Clock Model
  ClockModel cmodel = ClockModel::kMiniRafs;

  // measurements
  bool use_range = true;       // use GPS pseudorange measurement
  bool use_range_rate = true;  // use GPS pseudorange-rate measurement

  // Estimation
  int state_size = 8;          // Pos(3), vel(3), bias, drift [km, km/s, s, s/s]
  double pos_err = 1.0;        // Initial Position error [km]
  double vel_err = 1e-3;       // Initial Velocity error [km/s]
  double clk_bias_err = 1e-6;  // Initial Clock bias error [s]
  double clk_drift_err = 1e-9;  // Initial Clock drift error [s/s]
  double sigma_acc = 1e-12;     // Process noise Acceleration [km/s^2]  <-- tune
                                // this for optimal performance!

  // Debug mode
  bool plot_results = false;
  bool debug_jacobian = false;
  bool print_debug = false;
  bool no_meas = false;  // set to true to turn off measurements

  /**********************************************
   * Setup
   * ********************************************/

  // measurement type vector
  std::vector<GnssMeasurementType> meas_types;
  if (use_range) {
    meas_types.push_back(GnssMeasurementType::PR);
  }
  if (use_range_rate) {
    meas_types.push_back(GnssMeasurementType::PRR);
  }
  int meas_type_num = meas_types.size();

  // Orbit Dynamics
  auto dyn_earth_tb = std::make_shared<CartesianTwoBodyDynamics>(
      MU_EARTH);  // use 2d earth dynamics to propagate GPS constellation
  auto dyn_est = std::make_shared<NBodyDynamics>();   // Filter Dynamics
  auto dyn_true = std::make_shared<NBodyDynamics>();  // true dynamics

  auto earth = Body::Earth();
  auto moon_true = Body::Moon(moon_sph_true, moon_sph_true);
  auto moon_est = Body::Moon(moon_sph_est, moon_sph_est);

  dyn_true->AddBody(moon_true);
  dyn_true->SetCentralBody(moon_true);

  dyn_est->AddBody(moon_est);
  dyn_est->SetCentralBody(moon_est);

  if (add_earth) {
    dyn_true->AddBody(earth);
    dyn_est->AddBody(earth);
  }

  // clock dynamics
  auto dyn_clk = ClockDynamics(cmodel);

  // GPS constellation
  auto channel = std::make_shared<GnssChannel>();
  gps_const.SetChannel(channel);
  gps_const.SetDynamics(dyn_earth_tb);
  gps_const.LoadTleFile("gps");  // example gps file

  // Print Time
  std::string epoch_string = sp::TAItoStringUTC(epoch0, 3);
  std::cout << "Initial Epoch: " << epoch_string << std::endl;

  // Set dynamics integration time
  dyn_earth_tb->SetDt(dt);
  dyn_est->SetDt(dt);
  dyn_true->SetDt(dt);

  // Moon spacecraft
  ClassicalOE coe_moon({a, e, i, Omega, w, M}, CoordSystem::MI);
  auto cart_state_moon = std::make_shared<CartesianOrbitState>(
      ClassicalToCartesian(coe_moon, MU_MOON));

  Vector2 clock_vec{clk_bias, clk_drift};  // [s, s/s]
  ClockState clock_state(clock_vec);

  auto moon_sat = std::make_shared<Spacecraft>();
  auto receiver = std::make_shared<GnssReceiver>("moongpsr");

  moon_sat->AddDevice(receiver);
  moon_sat->SetDynamics(dyn_true);
  moon_sat->SetClock(clock_state);
  moon_sat->SetOrbitState(cart_state_moon);
  moon_sat->SetEpoch(epoch0);
  moon_sat->SetBodyId(NaifId::MOON);
  moon_sat->SetClockDynamics(dyn_clk);

  receiver->SetAgent(moon_sat);
  receiver->SetChannel(channel);
  channel->AddReceiver(receiver);

  // Initial covariance
  MatrixXd P0 =
      ConstructInitCovariance(pos_err, vel_err, clk_bias_err, clk_drift_err);

  // Joint state and dynamics
  JointState joint_state;
  joint_state.PushBackStateAndDynamics(cart_state_moon.get(), dyn_est.get());
  joint_state.PushBackStateAndDynamics(&clock_state, &dyn_clk);

  FilterDynamicsFunction joint_dynamics =
      joint_state.GetFilterDynamicsFunction();

  /*********************************************
   * Define Measurement function
   * *******************************************/
  FilterMeasurementFunction meas_func_pos_clk =
      [moon_sat, receiver, state_size, no_meas, meas_types](
          const VectorX x, MatrixXd& H, MatrixXd& R) -> VectorX {
    if (no_meas) {
      return VectorXd::Zero(0);
    }

    // Measurements
    std::string freq = "L1";
    double epoch = moon_sat->GetEpoch().val();
    auto measall =
        receiver->GetMeasurement(epoch);      // measurements of all frequencies
    auto meas = measall.ExtractSignal("L1");  // measurements of L1
    int sat_num =
        meas.GetTrackedSatelliteNum();       // number of tracked GPS satellites
    int mtot = sat_num * meas_types.size();  // total number of measurements

    // Predict measurements
    H = MatrixXd::Zero(mtot, x.size());
    R = MatrixXd::Zero(mtot, mtot);
    VectorX x_N = VectorX::Zero(mtot);  // a dummy variable for carrier phase
    CoordSystem coord_in = CoordSystem::MI;

    VectorX z = meas.GetPredictedGnssMeasurement(
        epoch, x.head(6), x.tail(2), x_N, H, meas_types,
        coord_in);  // Jacobian with autodiff

    // Get the Measurement Noise
    VectorXd noise_std_vec = meas.GetGnssNoiseStdVector(meas_types);
    R.diagonal().array() = noise_std_vec.array().square();

    return z;
  };

  /*********************************************
   * Define Process Noise function
   * *******************************************/
  FilterProcessNoiseFunction proc_noise_func = [cmodel, state_size, sigma_acc](
                                                   const VectorX& x,
                                                   real t_curr, real t_end) {
    int clock_index = 6;
    double dt = (t_end - t_curr).val();

    MatrixXd Q = MatrixXd::Zero(state_size, state_size);

    Matrix6d Q_rv = Matrix6d::Zero();
    for (int i = 0; i < 3; i++) {
      Q_rv(i, i) = pow(dt, 3) / 3.0 * pow(sigma_acc, 2);
      Q_rv(i + 3, i + 3) = dt * pow(sigma_acc, 2);
      Q_rv(i, i + 3) = pow(dt, 2) / 2.0 * pow(sigma_acc, 2);
      Q_rv(i + 3, i) = pow(dt, 2) / 2.0 * pow(sigma_acc, 2);
    }
    Matrix2d Q_clk = GetClockProcessNoise(cmodel, dt);

    Q.block(0, 0, 6, 6) = Q_rv;
    Q.block(6, 6, 2, 2) = Q_clk;

    return Q;
  };

  /*************************************
   * EKF Setup
   * ***********************************/
  EKF ekf;
  ekf.SetDynamicsFunction(joint_dynamics);
  ekf.SetMeasurementFunction(meas_func_pos_clk);
  ekf.SetProcessNoiseFunction(proc_noise_func);
  std::cout << "Initialized EKF" << std::endl;

  // Storage
  MatrixXd error_mat(4, time_step_num);
  VectorXd num_meas(time_step_num);

  // Initilization
  VectorX x_est = SampleMVN(joint_state.GetJointStateValue(), P0, 1, seed);
  ekf.Initialize(x_est, P0);
  VectorXd est_err = ComputeEstimationErrors(moon_sat, &ekf);
  error_mat.col(0) = est_err;

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

  /***********************************************
   * Main loop
   **********************************************/
  real t = t0;
  double epoch = epoch0;
  PrintProgressHeader();

  int time_index = 0;

  for (t = t0; t < tf; t += Dt) {
    time_index += 1;
    epoch += Dt;  // first propagate to the next epoch

    // Propagate True State
    moon_sat->Propagate(epoch);
    gps_const.Propagate(epoch);

    // Get True Measurement
    VectorX z_true;
    auto measall = receiver->GetMeasurement(epoch);
    auto meas = measall.ExtractSignal("L1");
    int num_sat = meas.GetTrackedSatelliteNum();
    num_meas(time_index) = num_sat;

    if (!no_meas) {
      bool with_noise = true;
      z_true = meas.GetGnssMeasurement(meas_types, with_noise, seed);
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

    // Compute Estimation
    est_err =
        ComputeEstimationErrors(moon_sat, &ekf);  // pos, vel, clkb, clkd error
    error_mat.col(time_index) = est_err;

    // Print progress
    if (fmod(t.val(), print_every) < 1e-3) {
      PrintProgress(t.val(), est_err(0), est_err(1), est_err(2));
    }

    // print measurement residuals
    if ((print_debug) && (!no_meas)) {
      PrintEKFDebugInfo(&ekf);
    }
  }

  // Print Statistics
  if (!no_meas) {
    PrintEstimationStatistics(num_meas, error_mat, 0.3);  // use last 30%

    // Write data ------------------------------------------
    std::cout << "Simulation finished, saving data..." << std::endl;
    writer.WriteData(*data_history);
    std::cout << "Data saved to " << output_path << std::endl;

    // Plots -----------------------------------------------
    if (plot_results) {
      Plot3DTrajectory(data_history, "true");
      Plot3DTrajectory(data_history, "est");

      // Plot estimate for all states
      PlotState(data_history, "true");
      PlotState(data_history, "est");
    }
  }

  return 0;
}
