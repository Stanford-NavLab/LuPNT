/**
 * @file filters.cpp
 * @author Stanford NAVLAB
 * @brief Implementation of Filters
 * @version 0.1
 * @date 2023-09-09
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/numerics/filters/ekf.h"

#include "lupnt/core/asset_factory.h"

namespace lupnt {

  KalmanFilter::KalmanFilter(Config& config) : Filter(config) {}

  EKF::EKF(Config& config) : KalmanFilter(config) {
    if (config["outlier_threshold"]) outlier_threshold_ = config["outlier_threshold"].as<double>();
  }

  void EKF::Predict(Real t, const State* u) {
    // Checks
    LUPNT_CHECK(x_.size() > 0, "(Predict) State vector is not set", "EKF");
    LUPNT_CHECK(P_.size() > 0, "(Predict) Covariance matrix is not set", "EKF");
    LUPNT_CHECK(f_dyn_, "(Predict) Dynamics function not set", "EKF");
    LUPNT_CHECK(f_proc_, "(Predict) Process noise function not set", "EKF");

    // Process noise
    Q_ = f_proc_(x_, t_, t);
    LUPNT_CHECK(!Q_.hasNaN(), "(Predict) Process noise covariance has NaN", "EKF");

    if (use_process_noise_mapping_) {
      Q_ = G_ * Q_ * G_.transpose();
      LUPNT_CHECK(!Q_.hasNaN(), "(Predict) Mapped process noise covariance has NaN", "EKF");
    }

    // Dynamics
    x_ = f_dyn_(x_, t_, t, u, &F_);
    LUPNT_CHECK(!x_.hasNaN(), "(Predict) State has NaN", "EKF");
    LUPNT_CHECK(F_.cols() == x_.rows() && F_.rows() == x_.rows(),
                "(Predict) State and STM have different sizes", "EKF");

    P_bar_ = F_ * P_ * F_.transpose();
    P_ = P_bar_ + Q_;
    LUPNT_CHECK(!P_.hasNaN(), "(Predict) Covariance has NaN", "EKF");
    // LUPNT_CHECK((P_.diagonal().array() >= 0).all(), "Covariance has negative diagonal", "EKF");

    t_ = t;

    // Store
    x_prior_ = x_;
    P_prior_ = P_;
  }

  void EKF::SetOutlierThreshold(double outlier_threshold) {
    LUPNT_CHECK(outlier_threshold >= 0, "(SetOutlierThreshold) Outlier threshold must be positive",
                "EKF");
    outlier_threshold_ = outlier_threshold;
  }

  void EKF::SetConsiderStateCount(int n_consider) {
    LUPNT_CHECK(n_consider >= 0, "(SetConsiderStateCount) n_consider must be non-negative", "EKF");
    n_consider_ = n_consider;
  }

  void EKF::Update(const VecX& z_true) {
    LUPNT_CHECK(!z_true.hasNaN(), "(Update) True measurement has NaN", "EKF");
    LUPNT_CHECK(f_meas_, "(Update) Measurement function not set", "EKF");

    // Store
    z_true_ = z_true;
    x_post_ = x_;
    P_post_ = P_;

    // Size
    int N_state = x_.size();
    int N_meas = z_true_.size();
    if (N_meas == 0) return;  // no measurement, nothing to update

    z_prior_ = f_meas_(x_, &H_, &R_);
    LUPNT_CHECK(!z_prior_.hasNaN(), "(Update) Predicted measurement has NaN", "EKF");

    S_ = R_ + H_ * P_ * H_.transpose();  // Measurement information
    LUPNT_CHECK(!S_.hasNaN(), "(Update) Innovation covariance has NaN", "EKF");

    dz_ = z_true_ - z_prior_;
    LUPNT_CHECK(!dz_.hasNaN(), "(Update) Measurement residual has NaN", "EKF");

    // Outlier rejection
    if (use_custom_fault_detection_ && f_fault_det_) {
      f_fault_det_(this);
    } else {
      RemoveOutliers();
    }
    if (dz_.size() == 0) return;  // All measurements are outliers

    // Update step
    K_ = P_ * H_.transpose() * S_.inverse();  // Kalman gain
    LUPNT_CHECK(!K_.hasNaN(), "(Update) Kalman gain has NaN", "EKF");

    if (n_consider_ > 0) {
      // Schmidt-Kalman: never correct the trailing consider states. The Joseph-form
      // covariance update below remains valid for this non-optimal gain, so the
      // consider states' (co)variance still updates consistently even though their
      // mean does not move.
      LUPNT_CHECK(n_consider_ <= N_state, "(Update) n_consider exceeds state size", "EKF");
      K_.bottomRows(n_consider_).setZero();
    }

    dx_ = K_ * dz_;
    Sigma_dx_ = K_ * S_ * K_.transpose();
    x_ = x_ + dx_;
    LUPNT_CHECK(!x_.hasNaN(), "(Update) State has NaN", "EKF");

    MatXd I = MatXd::Identity(N_state, N_state);
    MatXd G = I - K_ * H_;
    P_ = G * P_ * G.transpose() + K_ * R_ * K_.transpose();  // Joseph form
    LUPNT_CHECK(!P_.hasNaN(), "(Update) Covariance has NaN", "EKF");
    LUPNT_CHECK((P_.diagonal().array() >= 0).all(), "(Update) Covariance has negative diagonal",
                "EKF");

    // Covariance inflation
    // double lambda = 0.0;
    // P_ = (1 + lambda) * P_;
    // P_ = G * P_;

    // Store
    x_post_ = x_;
    P_post_ = P_;
  }

  void EKF::RemoveOutliers() {
    VecXd ratio = dz_.array().abs() / S_.diagonal().array().sqrt();
    VecXb is_outlier = ratio.array() > outlier_threshold_;

    int N_state = x_.size();
    int N_meas = dz_.size();
    int N_meas_new = N_meas - is_outlier.count();

    if (N_meas_new == N_meas) return;

    // Logger::Info(fmt::format("Removing outliers {}/{}", N_meas_new, N_meas), "EKF", t_);
    // if (N_meas_new == 0) {
    //   Logger::Info("All measurements are outliers", "EKF", t_);
    //   Logger::Info("dz: " + fmt::to_string(dz_.array().abs().transpose()), "EKF", t_);
    //   Logger::Info("S diag: " + fmt::to_string(S_.diagonal().array().sqrt().transpose()), "EKF",
    //                t_);
    // }

    VecXd dz_new(N_meas_new);
    MatXd H_new(N_meas_new, N_state);
    MatXd R_new(N_meas_new, N_meas_new);

    int n = 0;
    for (int i = 0; i < N_meas; i++) {
      if (is_outlier[i]) continue;
      dz_new(n) = dz_(i);
      H_new.row(n) = H_.row(i);

      int m = 0;
      for (int j = 0; j < N_meas; j++) {
        if (is_outlier[j]) continue;
        R_new(n, m) = R_(i, j);
        m++;
      }

      n++;
    }

    // Update
    H_ = H_new;
    R_ = R_new;
    dz_ = dz_new;
    S_ = H_ * P_prior_ * H_.transpose() + R_;
  }

  /******************************************************************************
   * Smoothing functions
   ******************************************************************************/

  void EKF::InitializeLogger(int max_tidx) {
    max_tidx_ = max_tidx;
    x_prior_log_.resize(max_tidx_);
    P_prior_log_.resize(max_tidx_);
    x_pos_log_.resize(max_tidx_);
    P_pos_log_.resize(max_tidx_);
    stm_log_.resize(max_tidx_);
    sm_tidx_ = max_tidx_ - 1;
    time_log_.resize(max_tidx_);
    x_sm_.resize(max_tidx_);
    P_sm_.resize(max_tidx_);
  }

  void EKF::InitializeSmootherState() {
    x_sm_[sm_tidx_] = x_post_;
    P_sm_[sm_tidx_] = P_post_;
    time_log_[sm_tidx_] = t_;
  }

  void EKF::LogFilterEstimate(int tidx) {
    x_prior_log_[tidx] = x_prior_;
    P_prior_log_[tidx] = P_prior_;
    x_pos_log_[tidx] = x_post_;
    P_pos_log_[tidx] = P_post_;
    stm_log_[tidx] = F_;  // (tidx-1 -> tidx)
    time_log_[tidx] = t_;
  }

  // Smoother step Obtain x_tidx|N. and P_tidx|N
  void EKF::UpdateSmoother(int tidx) {
    // Compute gain
    int k = tidx;
    MatXd P_kk = P_pos_log_[k];
    MatXd P_k1k = P_prior_log_[k + 1];
    MatXd stm_k = stm_log_[k + 1];
    VecXd x_kk = x_pos_log_[k];
    VecXd x_k1k = x_prior_log_[k + 1];

    MatXd A_k = P_kk * stm_k.transpose() * P_k1k.inverse();

    // Update smoothed state
    x_sm_[k] = x_kk + A_k * (x_sm_[k + 1] - x_k1k);
    // Update smoothed covariance
    P_sm_[k] = P_kk + A_k * (P_sm_[k + 1] - P_k1k) * A_k.transpose();
    sm_tidx_ = tidx;

    // Log time
    t_ = time_log_[k];
  }

  REGISTER_FACTORY_CLASS(Filter, EKF)

  // Define the GetRegistry function for this specialization (must come before explicit
  // instantiation)
  template <> std::unordered_map<std::string, AssetFactory<Filter, Config&>::Creator>&
  AssetFactory<Filter, Config&>::GetRegistry() {
    return Registry();
  }

  // Explicit template instantiation to ensure single registry across library boundaries
  template class AssetFactory<Filter, Config&>;

}  // namespace lupnt
