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
#include <lupnt/agents/agent.h>
#include <lupnt/agents/gnss_constellation.h>
#include <lupnt/agents/state_estimation_app.h>
#include <lupnt/core/file.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/measurements/gnss_channel.h>
#include <lupnt/measurements/gnss_measurement.h>
#include <lupnt/measurements/gnss_receiver.h>
#include <lupnt/measurements/transmission.h>
#include <lupnt/numerics/filters.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/clock.h>
#include <lupnt/physics/coord_converter.h>
#include <lupnt/physics/orbit_state.h>
#include <lupnt/physics/spice_interface.h>

// C++ includes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace lupnt;
namespace sp = SpiceInterface;

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
  double tf = t0 + 24.0 * SECS_PER_HOUR;  // 6 hours
  double dt = 1.0;                        // Integration time step [s]
  double Dt = 10.0;                       // Propagation time step [s]
  double print_every = 1.0 * SECS_PER_HOUR;
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
  Matrix6d P_rv = Matrix6d::Zero();
  P_rv.block(0, 0, 3, 3) = Matrix3d::Identity() * pow(pos_err, 2);
  P_rv.block(3, 3, 3, 3) = Matrix3d::Identity() * pow(vel_err, 2);

  Matrix2d P_clk = Matrix2d::Zero();
  P_clk(0, 0) = pow(clk_bias_err, 2);
  P_clk(1, 1) = pow(clk_drift_err, 2);

  MatrixXd P0 = blkdiag(P_rv, P_clk);

  // Joint state and dynamics
  JointState joint_state;
  joint_state.PushBackStateAndDynamics(cart_state_moon.get(),
                                       dyn_moon_tb.get());
  joint_state.PushBackStateAndDynamics(&clock_state, &dyn_clk);

  FilterDynamicsFunction joint_dynamics =
      joint_state.GetFilterDynamicsFunction();

  // Process noise
  FilterProcessNoiseFunction proc_noise_func =
      [state_size](const VectorX& x, real t_curr, real t_end) {
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

  VectorX x_est = SampleMVN(joint_state.GetJointStateValue(), P0, 1);
  ekf.Initialize(x_est, P0);

  // Output
  auto data_history = std::make_shared<DataHistory>();
  auto output_path = std::filesystem::current_path() / "output" / "ExampleEKF";
  FileWriter writer(output_path, true);

  // Main loop
  double t = t0;
  double epoch = epoch0;

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

    // Navigation
    data_history->AddData("z_true", t, ekf.z_true);
    data_history->AddData("z_pred", t, ekf.z_pred);
    data_history->AddData("CN0", t, meas.GetCN0());

    data_history->AddData("vis_earth", t, meas.GetEarthOccultation());
    data_history->AddData("vis_moon", t, meas.GetMoonOccultation());
    data_history->AddData("vis_antenna", t, meas.GetMoonOccultation());
    data_history->AddData("vis_ionos", t, meas.GetMoonOccultation());

    // Moon spacecraft
    auto state = moon_sat->GetCartesianGCRFStateAtEpoch(epoch);
    auto sate_mi = ConvertOrbitStateCoordSystem(state, epoch, CoordSystem::MI);
    auto state_gcrf =
        ConvertOrbitStateCoordSystem(state, epoch, CoordSystem::GCRF);
    data_history->AddData("rv_moon_mi", t, sate_mi->GetVector());
    data_history->AddData("rv_moon_gcrf", t, state_gcrf->GetVector());

    // Estimation
    data_history->AddData("rv", t, moon_sat->GetOrbitState()->GetVector());
    data_history->AddData("rv_pred", t, ekf.xbar.head(6));
    data_history->AddData("rv_est", t, ekf.x.head(6));

    data_history->AddData("clk", t, moon_sat->GetClockState().GetVector());
    data_history->AddData("clk_pred", t, ekf.xbar.tail(2));
    data_history->AddData("clk_est", t, ekf.x.tail(2));

    data_history->AddData("P_rv", t, ekf.P.diagonal().segment(0, 6));
    data_history->AddData("P_clk", t, ekf.P.diagonal().segment(6, 2));

    // GPS constellation
    for (int i = 0; i < gps_const.GetNumSatellites(); i++) {
      auto sate =
          gps_const.GetSatellite(i)->GetCartesianGCRFStateAtEpoch(epoch);
      auto sate_mi = ConvertOrbitStateCoordSystem(sate, epoch, CoordSystem::MI);
      auto state_gcrf =
          ConvertOrbitStateCoordSystem(sate, epoch, CoordSystem::GCRF);

      std::string name = "sat" + std::to_string(i);
      data_history->AddData(name + "_mi", t, sate->GetVector());
      data_history->AddData(name + "_gcrf", t, state_gcrf->GetVector());
    }

    // Bodies
    data_history->AddData(
        "earth_mi", t,
        CoordConverter::Convert(epoch, VectorX::Zero(6), CoordSystem::GCRF,
                                CoordSystem::MI));
    data_history->AddData(
        "moon_gcrf", t,
        CoordConverter::Convert(epoch, VectorX::Zero(6), CoordSystem::MI,
                                CoordSystem::GCRF));

    // Print progress
    if (fmod(t, print_every) < 1e-3) {
      std::cout << "t = " << std::fixed << std::setprecision(1) << t / 3600.0
                << "/" << tf / 3600.0 << " hours (" << (t / tf * 100) << "%)"
                << "  Number of Measurements: " << meas_size << std::endl
                << std::setprecision(4);
    }
  }

  // Write data
  std::cout << "Simulation finished, saving data..." << std::endl;
  writer.WriteData(*data_history);
  std::cout << "Data saved to " << output_path << std::endl;
  return 0;
}
