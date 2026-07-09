/**
 * @file ukf.cc
 * @author Stanford NAV LAB
 * @brief  Unscented Kalman Filter
 * @version 0.1
 * @date 2024-12-30
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "lupnt/numerics/filters/ukf.h"

namespace lupnt {

  void UKF::InitializeUkfParams() {
    lambda_ = alpha_ * alpha_ * (n_x_ + kappa_) - n_x_;
    n_sigma_ = 2 * n_x_ + 1;

    w_m_.resize(n_sigma_);
    w_c_.resize(n_sigma_);

    // Weights for the 0-th sigma point
    w_m_(0) = lambda_ / (n_x_ + lambda_);
    w_c_(0) = lambda_ / (n_x_ + lambda_) + (1.0 - alpha_ * alpha_ + beta_);

    // Weights for the remaining sigma points
    for (int i = 1; i < n_sigma_; i++) {
      w_m_(i) = 1.0 / (2.0 * (n_x_ + lambda_));
      w_c_(i) = 1.0 / (2.0 * (n_x_ + lambda_));
    }
  }

  MatX UKF::ComputeSigmaPoints(const State& state, const MatXd& cov) {
    // chol_cov is the Cholesky of (n_x_ + lambda_) * cov
    MatXd chol_cov = ((n_x_ + lambda_) * cov).llt().matrixL();

    // sigma_points is [n_x_ x n_sigma_]
    MatX sigma_points(n_x_, n_sigma_);

    // 1) central sigma point
    sigma_points.col(0) = state;

    // 2) +/- columns
    for (int i = 0; i < n_x_; i++) {
      sigma_points.col(i + 1) = state + chol_cov.col(i);
      sigma_points.col(i + 1 + n_x_) = state - chol_cov.col(i);
    }

    return sigma_points;
  }

  State UKF::UnscentedTransform(const MatX& sigma_points, MatXd& cov_out, const MatXd* Q) {
    // Compute the mean of the sigma points
    State state_mean = sigma_points * w_m_;

    // Compute the covariance of the sigma points
    cov_out = MatXd::Zero(n_x_, n_x_);
    for (int i = 0; i < n_sigma_; i++) {
      VecXd diff = (sigma_points.col(i) - state_mean).cast<double>();
      cov_out += w_c_(i) * diff * diff.transpose();
    }
    if (Q != nullptr) {
      LUPNT_CHECK(Q->rows() == n_x_ && Q->cols() == n_x_, "Process noise covariance of wrong size",
                  "UKF");
      cov_out += *Q;
    }  // Add process noise if required

    return state_mean;
  }

  void UKF::Predict(Real t, const State* u) {
    // Initialize UKF parameters if not done yet
    if (w_m_.size() == 0 || w_c_.size() == 0) {
      InitializeUkfParams();
    }

    // 1) Generate sigma points from current (x_, P_)
    MatX sigma_points = ComputeSigmaPoints(x_, P_);

    // 2) Propagate each sigma point through dynamics
    MatX sigma_points_pred(n_x_, n_sigma_);
    for (int i = 0; i < n_sigma_; i++) {
      // f_dyn_ signature: f_dyn_(x, t_curr, t_end, const State* u, MatXd* F)
      // we pass nullptr for the Jacobian and control
      sigma_points_pred.col(i) = f_dyn_(sigma_points.col(i), t_, t, u, nullptr);
    }

    // 3) Compute predicted mean/cov + add process noise
    MatXd Q_proc = f_proc_(x_, t_, t);  // Q for process noise

    MatXd cov_pred(n_x_, n_x_);
    State state_pred = UnscentedTransform(sigma_points_pred, cov_pred, &Q_proc);

    // 4) Store prior
    x_prior_ = state_pred;
    P_prior_ = cov_pred;

    // 5) Update filter's main x_, P_
    x_ = state_pred;
    P_ = cov_pred;

    // 6) Advance time
    t_ = t;
  }

  void UKF::Update(const VecX& z_obs) {
    // Use current x_, P_ from the predict step
    MatX sigma_points = ComputeSigmaPoints(x_, P_);

    // 1) Transform sigma points through measurement function
    //    also gather dimension of measurement
    MatXd h_dummy, r_dummy;
    VecX meas_dummy = f_meas_(x_, &h_dummy, &r_dummy);
    int n_z = meas_dummy.size();
    int n_x = x_.size();

    MatXd meas_sigma_points(n_z, n_sigma_);
    std::vector<MatXd> R_store(n_sigma_);

    for (int i = 0; i < n_sigma_; i++) {
      MatXd R_i(n_z, n_z);
      MatXd H_dum(n_z, n_x);
      // Does not need jacobian
      meas_sigma_points.col(i) = f_meas_(sigma_points.col(i), &H_dum, &R_i).cast<double>();
      R_store[i] = R_i;
    }

    // 2) Compute measurement mean/cov
    VecXd meas_mean = VecXd::Zero(n_z);
    S_ = MatXd::Zero(n_z, n_z);

    for (int i = 0; i < n_sigma_; i++) {
      meas_mean += w_m_(i) * meas_sigma_points.col(i);
    }
    for (int i = 0; i < n_sigma_; i++) {
      VecXd delta_meas = meas_sigma_points.col(i) - meas_mean;
      S_ += w_c_(i) * (delta_meas * delta_meas.transpose());
    }
    // Add measurement noise (assuming it's the same R for all sigma points)
    S_ += R_store[0];  // or a chosen R

    // 3) Compute cross-covariance
    MatXd cross_cov = MatXd::Zero(n_x_, n_z);
    for (int i = 0; i < n_sigma_; i++) {
      VecXd delta_state = (sigma_points.col(i) - x_).cast<double>();
      VecXd delta_meas = meas_sigma_points.col(i) - meas_mean;
      cross_cov += w_c_(i) * (delta_state * delta_meas.transpose());
    }

    // 4) Kalman Gain
    K_ = cross_cov * S_.inverse();

    // 5) Update
    dy_ = z_obs - meas_mean;
    dx_ = K_ * dy_;
    x_post_ = x_ + dx_;
    P_post_ = P_ - K_ * S_ * K_.transpose();

    // 6) Store result
    x_ = x_post_;
    P_ = P_post_;
    R_ = R_store[0];  // store the R used for the update
    z_true_ = z_obs;
    z_prior_ = meas_mean;
  }

}  // namespace lupnt
