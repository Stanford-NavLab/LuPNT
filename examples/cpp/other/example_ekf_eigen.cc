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
#include <lupnt/core/file.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/measurements/gnss_channel.h>
#include <lupnt/measurements/gnss_measurement.h>
#include <lupnt/measurements/gnss_receiver.h>
#include <lupnt/measurements/transmission.h>
#include <lupnt/numerics/filters.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/clock.h>
#include <lupnt/physics/frame_converter.h>
#include <lupnt/physics/orbit_state.h>
#include <lupnt/physics/spice_interface.h>

// C++ includes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace lupnt;

int main() {
  // Dynamics
  auto dyn_earth_tb = std::make_shared<CartesianTwoBodyDynamics>(GM_EARTH);
  auto dyn_moon_tb = std::make_shared<CartesianTwoBodyDynamics>(GM_MOON);
  auto dyn_moon_nb = std::make_shared<NBodyDynamics>();

  auto earth = Body::Earth();
  auto moon = Body::Moon(5, 5);
  dyn_moon_nb->AddBody(moon);
  dyn_moon_nb->AddBody(earth);
  dyn_moon_nb->SetPrimaryBody(moon);

  auto dyn_clk = ClockDynamics(ClockModel::kMicrosemiCsac);

  // GPS constellation
  auto channel = std::make_shared<GnssChannel>();
  auto gps_const = GnssConstellation();
  gps_const.SetChannel(channel);
  gps_const.SetDynamics(dyn_earth_tb);
  gps_const.LoadTleFile("gps");

  // Time
  // double epoch0 = String2TAI("2001/04/06 07:51:28.788 UTC").val();
  double epoch0 = gps_const.GetEpoch();
  std::string epoch_string = TAItoStringUTC(epoch0, 3);
  std::cout << "Initial Epoch: " << epoch_string << std::endl;

  double t0 = 0;
  double tf = t0 + 24.0 * SECS_HOUR;  // 6 hours
  double dt = 1.0;                    // Integration time step [s]
  double Dt = 10.0;                   // Propagation time step [s]
  double print_every = 1.0 * SECS_HOUR;
  double save_every = Dt;

  // Moon spacecraft
  Real a = 6541.4;
  Real e = 0.6;
  Real i = 65.5 * RAD;
  Real Omega = 0.0 * RAD;
  Real w = 90.0 * RAD;
  Real M = 0.0 * RAD;
  ClassicalOE coe_moon({a, e, i, Omega, w, M});
  coe_moon.SetCoordSystem(Frame::MOON_CI);

  auto cart_state_moon =
      std::make_shared<CartesianOrbitState>(Classical2Cart(coe_moon, GM_MOON));
  auto moon_sat = std::make_shared<Spacecraft>();
  auto receiver = std::make_shared<GnssReceiver>("moongpsr");

  moon_sat->AddDevice(receiver);
  moon_sat->SetDynamics(dyn_moon_nb);
  moon_sat->SetOrbitState(cart_state_moon);
  moon_sat->SetEpoch(epoch0);
  moon_sat->SetBodyId(NaifId::MOON);
  moon_sat->SetClockDynamics(dyn_clk);

  receiver->SetAgent(moon_sat);
  receiver->SetChannel(channel);
  channel->AddReceiver(receiver);

  // Estimation
  int n_state = 8;
  double pos_err = 1.0;            // Position error [km]
  double vel_err = 1e-3;           // Velocity error [km/s]
  double clk_bias_err = 1e-6;      // Clock bias error [s]
  double clk_drift_err = 1e-9;     // Clock drift error [s/s]
  double sigma_acc = 2e-6;         // Process noise [km/s^2]
  double sigma_range = 5e-3;       // Range measurement noise [km]
  double sigma_range_rate = 1e-6;  // Range rate measurement noise [km/s]

  Vec6 rv = moon_sat->GetOrbitState()->GetVec();
  Vec2 clk = Vec2::Zero();

  // OrbitState transition matrices
  Mat6d Phi_rv;
  Mat2d Phi_clk;

  // Initial covariance
  Mat6d P_rv = Mat6d::Zero();
  for (int i = 0; i < 3; i++) {
    P_rv(i, i) = pow(pos_err, 2);
    P_rv(i + 3, i + 3) = pow(vel_err, 2);
  }

  Mat2d P_clk = Mat2d::Zero();
  P_clk(0, 0) = pow(clk_bias_err, 2);
  P_clk(1, 1) = pow(clk_drift_err, 2);

  // Initial estimates
  auto zero6 = Vec6d::Zero();
  Vec6 rv_est = rv + SampleMVN(zero6, P_rv, 1);

  auto zero2 = Vec2d::Zero();
  Vec2 clk_est = clk + SampleMVN(zero2, P_clk, 1);

  // Process noise
  Mat6d Q_rv = Mat6d::Zero();
  for (int i = 0; i < 3; i++) {
    Q_rv(i, i) = pow(Dt, 3) / 3.0 * pow(sigma_acc, 2);
    Q_rv(i + 3, i + 3) = Dt * pow(sigma_acc, 2);
    Q_rv(i, i + 3) = pow(Dt, 2) / 2.0 * pow(sigma_acc, 2);
    Q_rv(i + 3, i) = pow(Dt, 2) / 2.0 * pow(sigma_acc, 2);
  }
  Mat2d Q_clk = GetClockProcessNoise(ClockModel::kMicrosemiCsac, Dt);

  // Output
  auto data_history = std::make_shared<DataHistory>();
  auto output_path = std::filesystem::current_path() / "output" / "ExampleEKF";
  FileWriter writer(output_path, true);

  // Main loop
  double t = t0;
  double epoch = epoch0;

  VecXd P(n_state, n_state);
  P.setZero();
  P.block(0, 0, 6, 6) = P_rv;
  P.block(6, 6, 2, 2) = P_clk;

  VecXd Q(n_state, n_state);
  Q.setZero();
  Q.block(0, 0, 6, 6) = Q_rv;
  Q.block(6, 6, 2, 2) = Q_clk;

  while (t < tf) {
    t += Dt;
    epoch += Dt;

    // Propagate
    moon_sat->Propagate(epoch);
    gps_const.Propagate(epoch);

    rv = moon_sat->GetOrbitState()->GetVec();
    clk = moon_sat->GetClockState().GetVec();

    // Predict
    auto rv_pred = rv_est;
    auto clk_pred = clk_est;
    dyn_moon_tb->PropagateWithStm(rv_pred, t - Dt, t, dt, Phi_rv);

    Phi_clk << 1.0, Dt, 0.0, 1.0;
    clk_pred = Phi_clk * clk_pred;

    // OrbitState and covariance
    VecX x(n_state);
    x << rv_pred, clk_pred;

    VecXd Phi(n_state, n_state);
    Phi.setZero();
    Phi.block(0, 0, 6, 6) = Phi_rv;
    Phi.block(6, 6, 2, 2) = Phi_clk;

    P = Phi * P * Phi.transpose() + Q;

    // Measurements
    std::string freq = "L1";
    auto measall = receiver->GetMeasurement(epoch);
    auto meas = measall.ExtractSignal("L1");
    int n_meas = meas.GetTrackedSatelliteNum();
    auto CN0 = meas.GetCN0();
    auto z_pr = meas.GetPseudorange();

    // Predict measurements
    VecXd H_pr(n_meas, n_state);
    VecX z_pr_pred =
        meas.GetPredictedPseudorange(epoch, rv_pred, clk_pred, H_pr);

    VecXd R_pr = VecXd::Zero(n_meas, n_meas);
    R_pr.diagonal().array() = pow(sigma_range, 2);

    // Update
    if (n_meas > 0) {
      // R and H matrices
      VecXd H(n_meas, n_state);
      H = H_pr;

      VecXd R(n_meas, n_meas);
      R = R_pr;

      VecXd I = VecXd::Identity(n_state, n_state);
      VecX y = z_pr - z_pr_pred;
      VecXd S = H * P * H.transpose() + R;
      VecXd K = P * H.transpose() * S.inverse();
      // A = C*inv(B)
      // auto decompC(C);  // decompose C with a suiting
      // decomposition VecXd A =
      // decompC.transpose().solve(B.transpose()).transpose();

      VecX dx = K * y;
      x = x + dx;
      P = (I - K * H) * P * (I - K * H).transpose() +
          K * R * K.transpose();  // Joseph form

      // Unpack state and covariance
      rv_est = x.segment(0, 6);
      clk_est = x.segment(6, 2);
      // Vec6 error = rv - rv_est;
      // std::cout << "error: " << error.transpose() << std::endl;

    } else {
      rv_est = rv_pred;
      clk_est = clk_pred;
    }

    // Navigation
    data_history->AddData("z_pr", t, z_pr);
    data_history->AddData("z_pr_pred", t, z_pr_pred);
    data_history->AddData("z_pr_st", t, z_pr_pred);
    data_history->AddData("CN0", t, meas.GetCN0());

    data_history->AddData("vis_earth", t, meas.GetEarthOccultation());
    data_history->AddData("vis_moon", t, meas.GetMoonOccultation());
    data_history->AddData("vis_antenna", t, meas.GetMoonOccultation());
    data_history->AddData("vis_ionos", t, meas.GetMoonOccultation());

    // Moon spacecraft
    auto state = moon_sat->GetCartesianGCRFStateAtEpoch(epoch);
    auto sate_mi = ConvertOrbitStateFrame(state, epoch, Frame::MOON_CI);
    auto state_gcrf = ConvertOrbitStateFrame(state, epoch, Frame::GCRF);
    data_history->AddData("rv_moon_mi", t, sate_mi->GetVec());
    data_history->AddData("rv_moon_gcrf", t, state_gcrf->GetVec());

    // Estimation
    data_history->AddData("rv", t, rv);
    data_history->AddData("rv_pred", t, rv_pred);
    data_history->AddData("rv_est", t, rv_est);

    data_history->AddData("clk", t, clk);
    data_history->AddData("clk_pred", t, clk_pred);
    data_history->AddData("clk_est", t, clk_est);

    data_history->AddData("P_rv", t, P.diagonal().segment(0, 6));
    data_history->AddData("P_clk", t, P.diagonal().segment(6, 2));

    // GPS constellation
    for (int i = 0; i < gps_const.GetNumSatellites(); i++) {
      auto sate =
          gps_const.GetSatellite(i)->GetCartesianGCRFStateAtEpoch(epoch);
      auto sate_mi = ConvertOrbitStateFrame(sate, epoch, Frame::MOON_CI);
      auto state_gcrf = ConvertOrbitStateFrame(sate, epoch, Frame::GCRF);

      std::string name = "sat" + std::to_string(i);
      data_history->AddData(name + "_mi", t, sate->GetVec());
      data_history->AddData(name + "_gcrf", t, state_gcrf->GetVec());
    }

    // Bodies
    Vec6 v6;
    v6.setZero();
    data_history->AddData("earth_mi", t,
                          ConvertFrame(epoch, v6, Frame::GCRF, Frame::MOON_CI));
    data_history->AddData("moon_gcrf", t,
                          ConvertFrame(epoch, v6, Frame::MOON_CI, Frame::GCRF));

    // Print progress
    if (fmod(t, print_every) < 1e-3) {
      std::cout << "t = " << std::fixed << std::setprecision(1) << t / 3600.0
                << "/" << tf / 3600.0 << " hours (" << (t / tf * 100) << "%)"
                << "  Number of Measurements: " << n_meas << std::endl
                << std::setprecision(4);
    }
  }

  // Write data
  std::cout << "Simulation finished, saving data..." << std::endl;
  writer.WriteData(*data_history);
  std::cout << "Data saved to " << output_path << std::endl;
  return 0;
}