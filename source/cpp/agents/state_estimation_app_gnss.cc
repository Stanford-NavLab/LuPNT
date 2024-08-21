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

#include "lupnt/agents/state_estimation_app_gnss.h"

namespace lupnt {

  void GnssStateEstimationApp::Setup() {
    Vec6 rv = agent->GetState()->GetVec();
    Vec2 clk = Vec2::Zero();

    // Initial covariance
    Mat6d P_rv = Mat6d::Zero();
    for (int i = 0; i < 3; i++) {
      P_rv(i, i) = pow(pos_err, 2);
      P_rv(i + 3, i + 3) = pow(vel_err, 2);
    }
    auto P_rv_pred_only = P_rv;

    Mat2d P_clk = Mat2d::Zero();
    P_clk(0, 0) = pow(clk_bias_err, 2);
    P_clk(1, 1) = pow(clk_drift_err, 2);
    auto P_clk_pred_only = P_clk;

    // Initial estimates
    auto zero6 = Vec6d::Zero();
    rv_est = rv + SampleMVN(zero6, P_rv, 1);

    auto zero2 = Vec2d::Zero();
    clk_est = clk + SampleMVN(zero2, P_clk, 1);

    rv_pred_only = rv_est;
    clk_pred_only = clk_est;

    // Process noise
    Q_rv = Mat6d::Zero();
    for (int i = 0; i < 3; i++) {
      Q_rv(i, i) = pow(Dt, 3) / 3.0 * pow(sigma_acc, 2);
      Q_rv(i + 3, i + 3) = Dt * pow(sigma_acc, 2);
      Q_rv(i, i + 3) = pow(Dt, 2) / 2.0 * pow(sigma_acc, 2);
      Q_rv(i + 3, i) = pow(Dt, 2) / 2.0 * pow(sigma_acc, 2);
    }
    Q_clk = ClockDynamics::TwoStateNoise(ClockModel::kMicrosemiCsac, Dt).cast<double>();

    P = VecXd::Zero(n_state, n_state);
    P.setZero();
    P.block(0, 0, 6, 6) = P_rv;
    P.block(6, 6, 2, 2) = P_clk;

    P_pred_only = VecXd::Zero(n_state, n_state);
    P_pred_only.setZero();
    P_pred_only.block(0, 0, 6, 6) = P_rv;
    P_pred_only.block(6, 6, 2, 2) = P_clk;

    Q = VecXd::Zero(n_state, n_state);
    Q.setZero();
    Q.block(0, 0, 6, 6) = Q_rv;
    Q.block(6, 6, 2, 2) = Q_clk;
  }

  void GnssStateEstimationApp::Step(double t) {
    epoch = epoch0 + t;

    MatXd Phi_6(6, 6);
    Vec6 rv_pred = dyn->Propagate(rv_est, t - Dt, t, &Phi_6);
    Phi_rv = Phi_6;
    rv_pred_only = dyn->Propagate(rv_pred_only, t - Dt, t, &Phi_6);
    Phi_rv_pred_only = Phi_6;

    Phi_clk << 1.0, Dt, 0.0, 1.0;
    Phi_clk_pred_only << 1.0, Dt, 0.0, 1.0;
    Vec2 clk_pred = Phi_clk * clk_est;
    clk_pred_only = Phi_clk * clk_pred_only;

    // State and covariance
    VecX x(n_state);
    x << rv_pred, clk_pred;

    VecXd Phi(n_state, n_state);
    Phi.setZero();
    Phi.block(0, 0, 6, 6) = Phi_rv;
    Phi.block(6, 6, 2, 2) = Phi_clk;

    VecXd Phi_pred_only(n_state, n_state);
    Phi_pred_only.setZero();
    Phi_pred_only.block(0, 0, 6, 6) = Phi_rv_pred_only;
    Phi_pred_only.block(6, 6, 2, 2) = Phi_clk_pred_only;

    P = Phi * P * Phi.transpose() + Q;
    P_pred_only = Phi_pred_only * P_pred_only * Phi_pred_only.transpose() + Q;

    // Measurements
    auto meas = receiver->GetMeasurement(epoch).ExtractSignal("L1");
    int n_meas = meas.GetTrackedSatelliteNum();
    auto CN0 = meas.GetCN0();
    auto z_pr = meas.GetPseudorange();

    // Predict measurements
    VecXd H_pr(n_meas, n_state);
    VecX z_pr_pred = meas.GetPredictedPseudorange(epoch, rv_pred, clk_pred, H_pr);

    VecXd R_pr = VecXd::Zero(n_meas, n_meas);
    R_pr.diagonal().array() = pow(sigma_range, 2);

    // Update
    if (n_meas > 0) {
      // R and H matrices
      VecXd H(n_meas, n_state);
      H = H_pr;

      VecXd R(n_meas, n_meas);
      R = R_pr;

      VecX y = z_pr - z_pr_pred;
      auto I = VecXd::Identity(n_state, n_state);
      auto S = H * P * H.transpose() + R;
      auto K = P * H.transpose() * S.inverse();
      // A = C*inv(B)
      // auto decompC(C);  // decompose C with a suiting
      // decomposition VecXd A =
      // decompC.transpose().solve(B.transpose()).transpose();

      VecX dx = K * y;
      x = x + dx;
      P = (I - K * H) * P * (I - K * H).transpose() + K * R * K.transpose();  // Joseph form

      // Unpack state and covariance
      rv_est = x.segment(0, 6);
      clk_est = x.segment(6, 2);
    } else {
      rv_est = rv_pred;
      clk_est = clk_pred;
    }
  }
}  // namespace lupnt
