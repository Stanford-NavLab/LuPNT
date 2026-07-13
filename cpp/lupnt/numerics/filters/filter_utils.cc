/**
 * @file filter_utils.cc
 * @author Stanford NAV Lab
 * @brief Utility functions for filters
 * @version 0.1
 * @date 2024-10-30
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "lupnt/numerics/filters/filter_utils.h"

#include "lupnt/dynamics/clock_dynamics.h"
#include "lupnt/numerics/filters/filter.h"

namespace lupnt {

  FilterMeasurementFunction GetFilterMeasurementFunction(const MeasurementFunction& f_meas) {
    FilterMeasurementFunction f_meas_func = [f_meas](const State& x, MatXd* H, MatXd* R) -> VecX {
      VecX y;
      VecX x_tmp = x.cast<double>();
      jacobian(f_meas, wrt(x_tmp), at(x_tmp, R), y, *H);
      return y;
    };
    return f_meas_func;
  }

  FilterDynamicsFunction GetFilterDynamicsFunction(const DynamicsFunction& f_dyn) {
    FilterDynamicsFunction f_dyn_func
        = [f_dyn](const State& x, Real t0, Real tf, const State* u, MatXd* F) -> State {
      VecX xf;
      VecX x_tmp = x.cast<double>();
      jacobian(f_dyn, wrt(x_tmp), at(x_tmp, t0, tf, u), xf, *F);
      return State(xf, x.GetNames(), x.GetUnits());
    };
    return f_dyn_func;
  }

  ProcessNoiseFunction GetProcessNoiseFunction(const ProcessNoiseFunction& f_proc) {
    // Identity pass-through: a ProcessNoiseFunction already has the filter
    // signature, so no adaptation is needed (mirrors the declared contract).
    return f_proc;
  }

  MatXd InitialCovariancePosVelClock(double sigma_r, double sigma_v, double sigma_b,
                                     double sigma_d) {
    Mat6d P_rv = Mat6d::Zero();
    P_rv.block(0, 0, 3, 3) = Mat3d::Identity() * pow(sigma_r, 2);
    P_rv.block(3, 3, 3, 3) = Mat3d::Identity() * pow(sigma_v, 2);

    Mat2d P_clk = Mat2d::Zero();
    P_clk(0, 0) = pow(sigma_b, 2);
    P_clk(1, 1) = pow(sigma_d, 2);

    MatXd P0 = BlockDiagonal(P_rv, P_clk);
    return P0;
  };

  Vec3d ProcessNoisePosVelCoeffs(Real dt) {
    LUPNT_CHECK(dt > 0.0, "dt must be set", "Asnc");
    Vec3d C_coeffs{pow(dt, 3) / 3.0, pow(dt, 2) / 2.0, dt};  // [C_11, C_21, C_22]
    return C_coeffs;
  }

  Mat6Xd ProcessNoisePosVelAccCoeffs(Real dt, const VecXd& beta) {
    LUPNT_CHECK(dt > 0.0, "dt must be set", "Admc");
    LUPNT_CHECK(beta.size() > 0, "beta must be set", "Admc");
    int n = beta.size();

    ArrXd b = beta.array();
    ArrXd b_inv = 1.0 / beta.array();
    ArrXd b_inv2 = b_inv.pow(2);
    ArrXd b_inv3 = b_inv.pow(3);
    ArrXd b_inv4 = b_inv.pow(4);
    ArrXd b_inv5 = b_inv.pow(5);

    ArrXd expbt = (-b * dt).exp();
    ArrXd expbt2 = (-2.0 * b * dt).exp();

    double dt2 = dt * dt;
    double dt3 = dt2 * dt;

    Mat6Xd C_coeffs(6, n);  // Coefficients diagonals

    // C_11
    C_coeffs.row(0) = (dt3 / 3.0) * b_inv2 - dt2 * b_inv3 + dt * b_inv4 * (1. - 2. * expbt)
                      + 0.5 * b_inv5 * (1. - expbt2);
    // C_21
    C_coeffs.row(1) = 0.5 * dt2 * b_inv2 - dt * b_inv3 * (1. - expbt) + b_inv4 * (1. - expbt)
                      - 0.5 * b_inv4 * (1. - expbt2);
    // C_22
    C_coeffs.row(2) = dt * b_inv2 - 2.0 * b_inv3 * (1. - expbt) + 0.5 * b_inv3 * (1. - expbt2);

    // C_31
    C_coeffs.row(3) = 0.5 * b_inv3 * (1 - expbt2) - dt * b_inv2 * expbt;
    // C_32
    C_coeffs.row(4) = 0.5 * b_inv2 * (1. + expbt2) - b_inv2 * expbt;
    // C_33
    C_coeffs.row(5) = 0.5 * b_inv * (1. - expbt2);
    return C_coeffs;
  }

  MatXd ProcessNoisePosVel(const MatXd& Q_a, Real dt) {
    const int n = Q_a.rows();
    Vec3d C_coeffs = ProcessNoisePosVelCoeffs(dt);
    MatXd Q_a_mat = (Q_a.rows() == Q_a.cols()) ? Q_a : Q_a.asDiagonal();

    MatXd Q_rv = MatXd::Zero(n * 2, n * 2);
    Q_rv.block(0, 0, n, n) = Q_a_mat * C_coeffs(0);
    Q_rv.block(0, n, n, n) = Q_a_mat * C_coeffs(1);
    Q_rv.block(n, 0, n, n) = Q_a_mat * C_coeffs(1);
    Q_rv.block(n, n, n, n) = Q_a_mat * C_coeffs(2);
    return Q_rv;
  }

  MatXd ProcessNoiseMappingMatrixPosVel(Real dt) {
    // Placeholder for future implementation
    MatXd M_rv = MatXd::Zero(6, 3);
    Vec3d C_coeffs = ProcessNoisePosVelCoeffs(dt);
    M_rv.block(0, 0, 3, 3) = C_coeffs(0) * Mat3d::Identity();
    M_rv.block(0, 3, 3, 3) = C_coeffs(1) * Mat3d::Identity();
    M_rv.block(3, 3, 3, 3) = C_coeffs(2) * Mat3d::Identity();
    return M_rv;
  }

  MatXd ProcessNoisePosVelAcc(const MatXd& Q_a, Real dt, const VecXd& beta) {
    const int n = Q_a.rows();
    Mat6Xd C_coeffs = ProcessNoisePosVelAccCoeffs(dt, beta);
    MatXd Q_a_mat = (Q_a.rows() == Q_a.cols()) ? Q_a : Q_a.asDiagonal();

    VecXd C11 = C_coeffs.row(0);
    VecXd C21 = C_coeffs.row(1);
    VecXd C22 = C_coeffs.row(2);
    VecXd C31 = C_coeffs.row(3);
    VecXd C32 = C_coeffs.row(4);
    VecXd C33 = C_coeffs.row(5);

    MatXd Q = MatXd::Zero(3 * n, 3 * n);
    Q.block(0, 0, n, n) = C11.asDiagonal() * Q_a_mat;
    Q.block(n, n, n, n) = C22.asDiagonal() * Q_a_mat;
    Q.block(2 * n, 2 * n, n, n) = C33.asDiagonal() * Q_a_mat;

    Q.block(0, n, n, n) = C21.asDiagonal() * Q_a_mat;
    Q.block(n, 0, n, n) = Q.block(0, n, n, n).transpose();

    Q.block(0, 2 * n, n, n) = C31.asDiagonal() * Q_a_mat;
    Q.block(2 * n, 0, n, n) = Q.block(0, 2 * n, n, n).transpose();
    Q.block(n, 2 * n, n, n) = C32.asDiagonal() * Q_a_mat;
    Q.block(2 * n, n, n, n) = Q.block(n, 2 * n, n, n).transpose();

    LUPNT_CHECK((Q.diagonal().array() >= 0).all(), "Q has negative diagonal", "AdaptiveModel");
    return Q;
  }

  MatXd ProcessNoiseClock(ClockModel clock_model, int n_clk, Real dt) {
    MatXd Q_clk(2, 2);

    if (n_clk == 2) {
      Q_clk.resize(2, 2);
      Q_clk = ClockDynamics::TwoStateNoise(clock_model, dt).cast<double>();
      return Q_clk;
    }

    if (n_clk == 3) {
      Q_clk.resize(3, 3);
      Q_clk = ClockDynamics::ThreeStateNoise(clock_model, dt).cast<double>();
      return Q_clk;
    }

    LUPNT_CHECK(false, "Invalid clock state size: must be 2 or 3", "ProcessNoiseClock");
    return Q_clk;
  }

  MatXd StateTransitionMatrixPosVel(double dt, int N_dim) {
    MatXd I = MatXd::Identity(N_dim, N_dim);
    MatXd Phi = MatXd::Zero(2 * N_dim, 2 * N_dim);
    Phi.block(0, 0, N_dim, N_dim) = I;
    Phi.block(N_dim, N_dim, N_dim, N_dim) = I;
    Phi.block(0, N_dim, N_dim, N_dim) = I * dt;
    return Phi;
  }

  MatXd StateTransitionMatrixPosVelAcc(double dt, const VecXd& beta) {
    const int n = beta.size();
    ArrXd b = beta.array();
    ArrXd b_inv = 1.0 / beta.array();
    ArrXd b_inv2 = b_inv.pow(2);
    ArrXd exp_bt = (-b * dt).exp();

    MatXd Phi = MatXd::Zero(3 * n, 3 * n);
    MatXd Phi_rv = StateTransitionMatrixPosVel(dt, n);
    MatXd Phi_a = exp_bt.matrix().asDiagonal();
    MatXd Phi_ra = (b_inv * dt + b_inv2 * (exp_bt - 1)).matrix().asDiagonal();
    MatXd Phi_va = (b_inv * (1 - exp_bt)).matrix().asDiagonal();

    Phi.block(0, 0, 2 * n, 2 * n) = Phi_rv;
    Phi.block(0, 2 * n, n, n) = Phi_ra;
    Phi.block(n, 2 * n, n, n) = Phi_va;
    Phi.block(2 * n, 2 * n, n, n) = Phi_a;
    return Phi;
  }

  MatXd ProcessNoisePosVelClock(ClockModel clock_model, int n_clk, double sigma_a, int n_sat,
                                Real dt) {
    MatXd Q = MatXd::Zero(n_clk * n_sat, n_clk * n_sat);
    MatXd Q_a = sigma_a * sigma_a * MatXd::Identity(3, 3);
    MatXd Q_rv = ProcessNoisePosVel(Q_a, dt);

    for (int k = 0; k < n_sat; k++) {
      Mat2d Q_clk = ClockDynamics::TwoStateNoise(clock_model, dt).cast<double>();
      Q.block(k * n_clk, k * n_clk, 6, 6) = Q_rv;
      Q.block(k * n_clk + 6, k * n_clk + 6, 2, 2) = Q_clk;
    }
    return Q;
  }

  VecXd ComputeEstimationErrorPVC(const Ptr<Agent>& sat, Filter* filter, int start_idx = 0) {
    auto x_est = filter->GetStatePost();
    auto x_true = sat->GetState();

    double x_pos_err = 1000 * (x_true.segment(0, 3) - x_est.segment(start_idx, 3)).norm().val();
    double x_vel_err = 1e6 * (x_true.segment(3, 3) - x_est.segment(start_idx + 3, 3)).norm().val();
    double x_clk_bias_err = 3e8 * abs((x_true(6) - x_est(start_idx + 6)).val());
    double x_clk_drift_err = 3e8 * abs((x_true(7) - x_est(start_idx + 7)).val());

    VecXd est_err(4);
    est_err << x_pos_err, x_vel_err, x_clk_bias_err, x_clk_drift_err;

    return est_err;
  }

  VecXd ComputeEstimationErrorPVC(const std::vector<Ptr<Agent>>& sats, Filter* filter) {
    VecXd est_err = VecXd::Zero(4 * sats.size());

    for (size_t i = 0; i < sats.size(); i++) {
      VecXd est_err_i = ComputeEstimationErrorPVC(sats[i], filter, i * 8);
      est_err.segment(i * 4, 4) = est_err_i;
    }

    return est_err;
  }

  VecX ConstructTrueStateVecFromSats(const std::vector<Ptr<Agent>>& sats) {
    VecX x_true = VecX::Zero(8 * sats.size());
    for (size_t i = 0; i < sats.size(); i++) {
      VecX x_sat = sats[i]->GetState();
      x_true.segment(i * 8, 8) = x_sat;
    }

    return x_true;
  }

}  // namespace lupnt
