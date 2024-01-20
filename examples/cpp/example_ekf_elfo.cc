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

#include "example_utils.cc"

using namespace lupnt;
namespace sp = SpiceInterface;

// Util Functions

int main() {
  /**********************************************
   * Simulation Parameters
   *********************************************/
  // Time
  double t0 = 0;
  double tf = t0 + 12.0 * SECS_PER_HOUR;  // 6 hours
  double dt = 1.0;                        // Integration time step [s]
  double Dt = 10.0;                       // Propagation time step [s]
  double print_every = 600;
  double save_every = Dt;

  // Simulation seed
  int seed = 5;
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

  // Dynamics Model
  int true_moon_sph =
      5;  // degree and order of spherical harmonics in true dynamics

  // Onboard Clock
  ClockModel cmodel = ClockModel::kMiniRafs;

  // measurements
  bool use_range = true;
  bool use_range_rate = true;

  // Estimation
  int state_size = 8;           // Pos, vel, bias, drift [km, km/s, s, s/s]
  double pos_err = 1.0;         // Position error [km]
  double vel_err = 1e-3;        // Velocity error [km/s]
  double clk_bias_err = 1e-6;   // Clock bias error [s]
  double clk_drift_err = 1e-9;  // Clock drift error [s/s]
  double sigma_acc =
      1e-9;  // Process noise Acceleration [km/s^2]  <-- tune this!

  // Debug mode
  bool debug_jacobian = false;
  bool print_debug = false;
  bool no_meas = false;  // set to true to turn off measurements
  bool use_nbody = false;

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

  // Dynamics
  auto dyn_earth_tb = std::make_shared<CartesianTwoBodyDynamics>(MU_EARTH);
  auto dyn_moon_tb = std::make_shared<CartesianTwoBodyDynamics>(
      MU_MOON);  // Todo: update this to more high fidelity dynamics
  auto dyn_moon_nb = std::make_shared<NBodyDynamics>();

  if (use_nbody) {
    auto earth = Body::Earth();
    auto moon = Body::Moon(true_moon_sph, true_moon_sph);
    dyn_moon_nb->AddBody(moon);
    dyn_moon_nb->AddBody(earth);
    dyn_moon_nb->SetCentralBody(moon);
  } else {
    dyn_moon_nb->SetCentralBody(Body::Moon());
    dyn_moon_nb->AddBody(Body::Moon());
  }

  auto dyn_clk = ClockDynamics(cmodel);

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
  joint_state.PushBackStateAndDynamics(cart_state_moon.get(),
                                       dyn_moon_tb.get());
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
  int time_step_num = int((tf - t0) / Dt) + 1;
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
    auto CN0 = meas.GetCN0();
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
  PrintEstimationStatistics(num_meas, error_mat, 0.3);  // use last 30%

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
