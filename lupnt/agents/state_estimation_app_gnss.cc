/**
 * @file GnssStateEstimationApp.cpp
 * @author Stanford NAVLAB
 * @brief State Estimation Application Using Gnss Measurement
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "state_estimation_app_gnss.h"

namespace lupnt {

void GnssStateEstimationApp::Setup() {
  ad::Vector6real rv = agent->GetOrbitState()->GetVector();
  ad::Vector2real clk = ad::Vector2real::Zero();

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
  rv_est = rv + SampleMVN(zero6, P_rv, 1);

  auto zero2 = Eigen::Vector2d::Zero();
  clk_est = clk + SampleMVN(zero2, P_clk, 1);

  rv_pred_only = rv_est;
  clk_pred_only = clk_est;

  // Process noise
  Q_rv = Eigen::Matrix6d::Zero();
  for (int i = 0; i < 3; i++) {
    Q_rv(i, i) = pow(Dt, 3) / 3.0 * pow(sigma_acc, 2);
    Q_rv(i + 3, i + 3) = Dt * pow(sigma_acc, 2);
    Q_rv(i, i + 3) = pow(Dt, 2) / 2.0 * pow(sigma_acc, 2);
    Q_rv(i + 3, i) = pow(Dt, 2) / 2.0 * pow(sigma_acc, 2);
  }
  Q_clk = GetClockProcessNoise(ClockModel::kMicrosemiCsac, Dt);

  P = Eigen::MatrixXd::Zero(n_state, n_state);
  P.setZero();
  P.block(0, 0, 6, 6) = P_rv;
  P.block(6, 6, 2, 2) = P_clk;

  P_pred_only = Eigen::MatrixXd::Zero(n_state, n_state);
  P_pred_only.setZero();
  P_pred_only.block(0, 0, 6, 6) = P_rv;
  P_pred_only.block(6, 6, 2, 2) = P_clk;

  Q = Eigen::MatrixXd::Zero(n_state, n_state);
  Q.setZero();
  Q.block(0, 0, 6, 6) = Q_rv;
  Q.block(6, 6, 2, 2) = Q_clk;
}

void GnssStateEstimationApp::Step(double t) {
  epoch = epoch0 + t;

  auto rv_pred = rv_est;
  auto clk_pred = clk_est;
  dyn->PropagateWithStm(rv_pred, t - Dt, t, dt, Phi_rv);
  dyn->PropagateWithStm(rv_pred_only, t - Dt, t, dt, Phi_rv_pred_only);
  Phi_clk << 1.0, Dt, 0.0, 1.0;
  Phi_clk_pred_only << 1.0, Dt, 0.0, 1.0;
  clk_pred = Phi_clk * clk_pred;
  clk_pred_only = Phi_clk * clk_pred_only;

  // State and covariance
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
  auto meas = receiver->GetMeasurement(epoch).ExtractSignal("L1");
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

    ad::VectorXreal y = z_pr - z_pr_pred;
    auto I = Eigen::MatrixXd::Identity(n_state, n_state);
    auto S = H * P * H.transpose() + R;
    auto K = P * H.transpose() * S.inverse();
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
  } else {
    rv_est = rv_pred;
    clk_est = clk_pred;
  }

  // Save data
  data_history->AddData("z_pr", t, z_pr);
  data_history->AddData("z_pr_pred", t, z_pr_pred);
  data_history->AddData("z_pr_st", t, z_pr_pred);
  data_history->AddData("CN0", t, meas.GetCN0());

  data_history->AddData("vis_earth", t, meas.GetEarthOccultation());
  data_history->AddData("vis_moon", t, meas.GetMoonOccultation());
  data_history->AddData("vis_antenna", t, meas.GetMoonOccultation());
  data_history->AddData("vis_ionos", t, meas.GetMoonOccultation());

  data_history->AddData("rv", t, agent->GetOrbitState()->GetVector());
  data_history->AddData("rv_pred", t, rv_pred);
  data_history->AddData("rv_pred_only", t, rv_pred_only);
  data_history->AddData("rv_est", t, rv_est);

  data_history->AddData("clk", t, agent->GetClock());
  data_history->AddData("clk_pred", t, clk_pred);
  data_history->AddData("clk_est", t, clk_est);
  data_history->AddData("clk_pred_only", t, clk_pred_only);

  data_history->AddData("P_rv", t, P.diagonal().segment(0, 6));
  data_history->AddData("P_rv_pred_only", t,
                        P_pred_only.diagonal().segment(0, 6));
  data_history->AddData("P_clk", t, P.diagonal().segment(6, 2));
  data_history->AddData("P_clk_pred_only", t,
                        P_pred_only.diagonal().segment(6, 2));
}
}  // namespace lupnt