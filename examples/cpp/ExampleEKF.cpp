// lupnt includes
#include <lupnt/agents/Agent.h>
#include <lupnt/agents/GNSSConstellation.h>
#include <lupnt/core/File.h>
#include <lupnt/dynamics/Dynamics.h>
#include <lupnt/measurements/GNSSChannel.h>
#include <lupnt/measurements/GNSSMeasurement.h>
#include <lupnt/measurements/GNSSReceiver.h>
#include <lupnt/measurements/Transmission.h>
#include <lupnt/numerics/Filters.h>
#include <lupnt/numerics/MathUtils.h>
#include <lupnt/numerics/eigenmvn.h>
#include <lupnt/physics/Clock.h>
#include <lupnt/physics/CoordConverter.h>
#include <lupnt/physics/OrbitState.h>
#include <lupnt/physics/SpiceInterface.h>

// Autodiff includes
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/QR>

// C++ includes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace LPT;
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
  auto channel = std::make_shared<GNSSChannel>();
  auto gps_const = GNSSConstellation();
  gps_const.SetChannel(channel);
  gps_const.SetDynamics(dyn_earth_tb);
  gps_const.LoadTleFile("gps");

  // Time
  // double epoch0 = sp::StringToTAI("2001/04/06 07:51:28.788 UTC").val();
  double epoch0 = gps_const.GetEpoch();
  std::string epoch_string = sp::TAItoStringUTC(epoch0, 3);
  std::cout << "Initial Epoch: " << epoch_string << std::endl;

  double t0 = 0;
  double tf = t0 + 6.0 * SECS_PER_HOUR;  // 6 hours
  double dt = 1.0;                       // Integration time step [s]
  double Dt = 10.0;                      // Propagation time step [s]
  double print_every = 1.0 * SECS_PER_HOUR;
  double save_every = Dt;

  // Moon spacecraft
  ad::real a = 6541.4;
  ad::real e = 0.6;
  ad::real i = 65.5 * RAD_PER_DEG;
  ad::real Omega = 90 * RAD_PER_DEG;
  ad::real w = 0.0 * RAD_PER_DEG;
  ad::real M = 0.0 * RAD_PER_DEG;
  ClassicalOE coe_moon(a, e, i, Omega, w, M);
  coe_moon.SetCoordSystem(CoordSystem::MI);

  auto cart_state_moon =
      std::make_shared<CartesianOrbitState>(CoeToCart(coe_moon, MU_MOON));
  auto moon_sat = std::make_shared<Spacecraft>();
  auto receiver = std::make_shared<GNSSReceiver>("moongpsr");

  moon_sat->AddDevice(receiver);
  moon_sat->SetDynamics(dyn_moon_nb);
  moon_sat->SetOrbitState(cart_state_moon);
  moon_sat->SetEpoch(epoch0);
  moon_sat->SetBodyId(BodyId::MOON);
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

  ad::Vector6real rv = moon_sat->GetOrbitState()->GetVector();
  ad::Vector2real clk = ad::Vector2real::Zero();

  // OrbitState transition matrices
  Eigen::Matrix6d Phi_rv;
  Eigen::Matrix6d Phi_rv_pred_only;
  Eigen::Matrix2d Phi_clk;
  Eigen::Matrix2d Phi_clk_pred_only;

  // Initial covariance
  Eigen::Matrix6d P_rv = Eigen::Matrix6d::Zero();
  for (int i = 0; i < 3; i++) {
    P_rv(i, i) = pow(pos_err, 2);
    P_rv(i + 3, i + 3) = pow(vel_err, 2);
  }
  auto P_rv_pred_only = P_rv;

  Eigen::Matrix2d P_clk = Eigen::Matrix2d::Zero();
  P_clk(0, 0) = pow(clk_bias_err, 2);
  P_clk(1, 1) = pow(clk_drift_err, 2);
  auto P_clk_pred_only = P_clk;

  // Initial estimates
  auto zero6 = Eigen::Vector6d::Zero();
  Eigen::EigenMultivariateNormal<double> mvrnd_rv(zero6, P_rv);
  ad::Vector6real rv_est = rv + mvrnd_rv.samples(1);

  auto zero2 = Eigen::Vector2d::Zero();
  Eigen::EigenMultivariateNormal<double> mvrnd_clk(zero2, P_clk);
  ad::Vector2real clk_est = clk;  // + mvrnd_clk.samples(1);

  auto rv_pred_only = rv_est;
  auto clk_pred_only = clk_est;

  // Process noise
  Eigen::Matrix6d Q_rv = Eigen::Matrix6d::Zero();
  for (int i = 0; i < 3; i++) {
    Q_rv(i, i) = pow(Dt, 3) / 3.0 * pow(sigma_acc, 2);
    Q_rv(i + 3, i + 3) = Dt * pow(sigma_acc, 2);
    Q_rv(i, i + 3) = pow(Dt, 2) / 2.0 * pow(sigma_acc, 2);
    Q_rv(i + 3, i) = pow(Dt, 2) / 2.0 * pow(sigma_acc, 2);
  }
  Eigen::Matrix2d Q_clk = GetClockProcessNoise(ClockModel::kMicrosemiCsac, Dt);

  // Output
  auto data_history = std::make_shared<DataHistory>();
  auto output_path = BASEPATH / "output" / "ExampleEKF";
  FileWriter writer(output_path, true);

  // Main loop
  double t = t0;
  double epoch = epoch0;

  Eigen::MatrixXd P(n_state, n_state);
  P.setZero();
  P.block(0, 0, 6, 6) = P_rv;
  P.block(6, 6, 2, 2) = P_clk;

  Eigen::MatrixXd P_pred_only(n_state, n_state);
  P_pred_only.setZero();
  P_pred_only.block(0, 0, 6, 6) = P_rv;
  P_pred_only.block(6, 6, 2, 2) = P_clk;

  Eigen::MatrixXd Q(n_state, n_state);
  Q.setZero();
  Q.block(0, 0, 6, 6) = Q_rv;
  Q.block(6, 6, 2, 2) = Q_clk;

  while (t < tf) {
    t += Dt;
    epoch += Dt;

    // Propagate
    moon_sat->Propagate(epoch);
    gps_const.Propagate(epoch);

    rv = moon_sat->GetOrbitState()->GetVector();
    clk = moon_sat->GetClock();

    // Predict
    auto rv_pred = rv_est;
    auto clk_pred = clk_est;
    dyn_moon_tb->PropagateWithStm(rv_pred, t - Dt, t, dt, Phi_rv);
    dyn_moon_tb->PropagateWithStm(rv_pred_only, t - Dt, t, dt,
                                  Phi_rv_pred_only);

    Phi_clk << 1.0, Dt, 0.0, 1.0;
    Phi_clk_pred_only << 1.0, Dt, 0.0, 1.0;
    clk_pred = Phi_clk * clk_pred;
    clk_pred_only = Phi_clk * clk_pred_only;

    // OrbitState and covariance
    ad::VectorXreal x(n_state);
    x << rv_pred, clk_pred;

    Eigen::MatrixXd Phi(n_state, n_state);
    Phi.setZero();
    Phi.block(0, 0, 6, 6) = Phi_rv;
    Phi.block(6, 6, 2, 2) = Phi_clk;

    Eigen::MatrixXd Phi_pred_only(n_state, n_state);
    Phi_pred_only.setZero();
    Phi_pred_only.block(0, 0, 6, 6) = Phi_rv_pred_only;
    Phi_pred_only.block(6, 6, 2, 2) = Phi_clk_pred_only;

    P = Phi * P * Phi.transpose() + Q;
    P_pred_only = Phi_pred_only * P_pred_only * Phi_pred_only.transpose() + Q;

    // Measurements
    std::string freq = "L1";
    auto measall = receiver->GetMeasurement(epoch);
    auto meas = measall.ExtractSignal("L1");
    int n_meas = meas.GetNumMeasurements();
    auto CN0 = meas.GetCN0();
    auto z_pr = meas.GetPseudorange();

    // Predict measurements
    Eigen::MatrixXd H_pr(n_meas, n_state);
    ad::VectorXreal z_pr_pred =
        meas.GetPseudorange(epoch, rv_pred, clk_pred, H_pr);

    Eigen::MatrixXd R_pr = Eigen::MatrixXd::Zero(n_meas, n_meas);
    R_pr.diagonal().array() = pow(sigma_range, 2);

    // Update
    if (n_meas > 0) {
      // R and H matrices
      Eigen::MatrixXd H(n_meas, n_state);
      H = H_pr;

      Eigen::MatrixXd R(n_meas, n_meas);
      R = R_pr;

      Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n_state, n_state);
      ad::VectorXreal y = z_pr - z_pr_pred;
      Eigen::MatrixXd S = H * P * H.transpose() + R;
      Eigen::MatrixXd K = P * H.transpose() * S.inverse();
      // A = C*inv(B)
      // auto decompC(C);  // decompose C with a suiting
      // decomposition Eigen::MatrixXd A =
      // decompC.transpose().solve(B.transpose()).transpose();

      ad::VectorXreal dx = K * y;
      x = x + dx;
      P = (I - K * H) * P * (I - K * H).transpose() +
          K * R * K.transpose();  // Joseph form

      // Unpack state and covariance
      rv_est = x.segment(0, 6);
      clk_est = x.segment(6, 2);
      // ad::Vector6real error = rv - rv_est;
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
    auto sate_mi = ConvertOrbitStateCoordSystem(state, epoch, CoordSystem::MI);
    auto state_gcrf =
        ConvertOrbitStateCoordSystem(state, epoch, CoordSystem::GCRF);
    data_history->AddData("rv_moon_mi", t, sate_mi->GetVector());
    data_history->AddData("rv_moon_gcrf", t, state_gcrf->GetVector());

    // Estimation
    data_history->AddData("rv", t, rv);
    data_history->AddData("rv_pred", t, rv_pred);
    data_history->AddData("rv_pred_only", t, rv_pred_only);
    data_history->AddData("rv_est", t, rv_est);

    data_history->AddData("clk", t, clk);
    data_history->AddData("clk_pred", t, clk_pred);
    data_history->AddData("clk_est", t, clk_est);
    data_history->AddData("clk_pred_only", t, clk_pred_only);

    data_history->AddData("P_rv", t, P.diagonal().segment(0, 6));
    data_history->AddData("P_rv_pred_only", t,
                          P_pred_only.diagonal().segment(0, 6));
    data_history->AddData("P_clk", t, P.diagonal().segment(6, 2));
    data_history->AddData("P_clk_pred_only", t,
                          P_pred_only.diagonal().segment(6, 2));

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
        CoordConverter::Convert(ad::VectorXreal::Zero(6), epoch,
                                CoordSystem::GCRF, CoordSystem::MI));
    data_history->AddData(
        "moon_gcrf", t,
        CoordConverter::Convert(ad::VectorXreal::Zero(6), epoch,
                                CoordSystem::MI, CoordSystem::GCRF));

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
