/**
 * @file filter.h
 * @author Stanford NAV LAB
 * @brief List of Filters
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <Eigen/QR>

#include "lupnt/numerics/filters/filter.h"

namespace lupnt {

  /// @brief Base class for linear/linearized Kalman-type filters (EKF, UDUEKF), adding
  /// Kalman-gain bookkeeping, fault detection, and a fixed-interval smoother on top of
  /// `Filter`.
  ///
  /// Stores the per-step Kalman quantities (`K_`, `S_`, `dz_`, `dx_`, `H_`, `F_`, ...)
  /// produced by `Predict`/`Update` so they can be inspected (e.g. by `filter_print.h`
  /// helpers or estimation-error utilities in `filter_utils.h`) and provides hooks
  /// (`InitializeLogger`/`LogFilterEstimate`/`UpdateSmoother`) for a Rauch-Tung-Striebel
  /// (RTS) smoother pass implemented by derived classes (e.g. `EKF`).
  class KalmanFilter : public Filter {
  protected:
    MatXd Q_;  // Process noise cov
    MatXd R_;  // Measurement noise cov

    MatXd S_;         // Innovation cov
    MatXd K_;         // Kalman gain
    MatXd G_;         // Process noise mapping matrix
    MatXd Sigma_dx_;  // State update covariance
    VecXd dz_;        // Measurement residual
    VecXd dx_;        // State update
    VecXd z_true_;    // Observed measurement
    VecXd z_prior_;   // Predicted measurement

    MatXd F_;  // State transition jacobian
    MatXd H_;  // Measurement jacobian

    bool use_process_noise_mapping_ = false;
    bool use_custom_fault_detection_ = false;

    // Smoothing
    int max_tidx_;  // Maximum time index
    int sm_tidx_;   // Current time index
    std::vector<double> time_log_;
    std::vector<VecXd> x_sm_;
    std::vector<MatXd> P_sm_;

    std::function<void(KalmanFilter*)> f_fault_det_;

  public:
    KalmanFilter() = default;

    /// @brief Construct a Kalman filter from a YAML config node (forwarded to `Filter`).
    KalmanFilter(Config& config);
    virtual ~KalmanFilter() = default;

    /// @brief Get the last state correction `dx = K * dz` applied by `Update`, size [n_x].
    VecXd GetStateCorrection() { return dx_; }
    /// @brief Get the observed measurement vector from the last `Update` call, size [n_z].
    VecXd GetTrueMeasurement() { return z_true_; }
    /// @brief Get the predicted measurement (from `f_meas_`) from the last `Update` call,
    /// size [n_z].
    VecXd GetPredictedMeasurement() { return z_prior_; }
    /// @brief Get the measurement residual `dz = z_true - z_prior` from the last `Update`
    /// call (after outlier removal), size [n_z].
    VecXd GetMeasurementResidual() { return dz_; }

    /// @brief Get the process noise covariance `Q_` used in the last `Predict` call,
    /// size [n_x x n_x].
    MatXd GetProcessNoise() { return Q_; }
    /// @brief Get the innovation covariance `S = H P H^T + R` from the last `Update` call,
    /// size [n_z x n_z].
    MatXd GetInnovationCov() { return S_; }
    /// @brief Get the measurement noise covariance `R_` from the last `Update` call,
    /// size [n_z x n_z].
    MatXd GetMeasurementCov() { return R_; }
    /// @brief Get the Kalman gain `K = P H^T S^-1` from the last `Update` call,
    /// size [n_x x n_z].
    MatXd GetKalmanGain() { return K_; }
    /// @brief Get the process noise mapping matrix `G_` (maps a reduced-dimension process
    /// noise onto the full state, used when `use_process_noise_mapping_` is enabled via
    /// `SetProcessNoiseMappingMatrix`).
    MatXd GetProcessNoiseMappingMatrix() { return G_; }
    /// @brief Same as `GetMeasurementCov`; get the measurement noise covariance `R_`.
    MatXd GetMeasurementNoiseCov() { return R_; }
    /// @brief Get the state-correction covariance `Sigma_dx = K S K^T` from the last
    /// `Update` call, size [n_x x n_x].
    MatXd GetStateCorrectionCov() { return Sigma_dx_; }
    /// @brief Get the state-transition (Jacobian) matrix `F_` from the last `Predict` call,
    /// size [n_x x n_x].
    MatXd GetStateJacobian() { return F_; }
    /// @brief Get the measurement Jacobian `H_` from the last `Update` call,
    /// size [n_z x n_x].
    MatXd GetMeasurementJacobian() { return H_; }
    /// @brief Get the smoothed state estimate at time index `tidx`, produced by
    /// `UpdateSmoother`.
    VecXd GetSmoothedState(int tidx) { return x_sm_[tidx]; }
    /// @brief Get the smoothed state covariance at time index `tidx`, produced by
    /// `UpdateSmoother`.
    MatXd GetSmoothedCovariance(int tidx) { return P_sm_[tidx]; }

    // Setters for fault detection
    /// @brief Override the measurement Jacobian `H_` (used by a custom fault-detection
    /// function to drop/reweight measurements before `Update`'s gain computation).
    void SetMeasurementJacobian(const MatXd& H) { H_ = H; }
    /// @brief Override the measurement noise covariance `R_` (used by a custom
    /// fault-detection function).
    void SetMeasurementNoiseCov(const MatXd& R) { R_ = R; }
    /// @brief Override the innovation covariance `S_` (used by a custom fault-detection
    /// function).
    void SetInnovationCov(const MatXd& S) { S_ = S; }
    /// @brief Override the measurement residual `dz_` (used by a custom fault-detection
    /// function, e.g. to remove outlier components before the gain computation).
    void SetMeasurementResidual(const VecXd& dz) { dz_ = dz; }

    /// @brief Get the current measurement dimension (number of rows of `dz_`).
    int GetMeasurementSize() { return dz_.rows(); }

    /// @brief Set the process noise mapping matrix `G_` and enable process-noise mapping,
    /// so that `Predict` maps `Q_` via `Q_ = G_ * Q_ * G_^T` before adding it to the
    /// propagated covariance.
    ///
    /// @param G  Process noise mapping matrix, size [n_x x n_q]
    void SetProcessNoiseMappingMatrix(const MatXd& G) {
      G_ = G;
      use_process_noise_mapping_ = true;
    }

    /// @brief Register a custom fault-detection/outlier-rejection callback, invoked by
    /// `Update` (in place of the default `RemoveOutliers`) after the measurement residual
    /// `dz_` and innovation covariance `S_` have been computed but before the Kalman gain
    /// is applied.
    ///
    /// The callback typically calls `SetMeasurementResidual`/`SetMeasurementJacobian`/
    /// `SetMeasurementNoiseCov`/`SetInnovationCov` to remove or reweight outlier components.
    void SetFaultDetectionFunction(std::function<void(KalmanFilter*)> f_fault_det) {
      f_fault_det_ = f_fault_det;
      use_custom_fault_detection_ = true;
    }

    /// @brief Allocate per-time-step smoother history buffers for `max_tidx_` epochs.
    ///
    /// Called once before a forward filtering + backward smoothing run to size the
    /// `x_prior_log_`/`P_prior_log_`/`x_pos_log_`/`P_pos_log_`/`stm_log_`/`x_sm_`/`P_sm_`
    /// buffers used by `LogFilterEstimate` and `UpdateSmoother`.
    ///
    /// @param max_tidx  Number of time steps to allocate for
    virtual void InitializeLogger(int max_tidx) = 0;  // Initialize logger for smoother

    /// @brief Seed the backward smoother pass with the final posterior state/covariance
    /// (`x_sm_[N-1] = x_post_`, `P_sm_[N-1] = P_post_`).
    ///
    /// Called once, after the forward filtering pass completes, before the first call to
    /// `UpdateSmoother`.
    virtual void InitializeSmootherState() = 0;  // Initialize smoother state

    /// @brief Record the prior/posterior state, covariance, and state-transition matrix at
    /// time index `tidx` for later use by `UpdateSmoother`.
    ///
    /// Called once per epoch during the forward filtering pass, after each `Predict`/
    /// `Update` pair.
    ///
    /// @param tidx  Time index at which to record the current filter state
    virtual void LogFilterEstimate(int tidx) = 0;  // Log filter estimate at time index

    /// @brief Perform one backward Rauch-Tung-Striebel (RTS) smoother step, computing the
    /// smoothed state/covariance at time index `tidx` from the smoothed values at `tidx+1`
    /// and the logged filter history.
    ///
    /// Called once per epoch, in decreasing time-index order, during the backward smoothing
    /// pass after the forward filtering pass and `InitializeSmootherState` have completed.
    ///
    /// @param tidx  Time index to smooth (must be `< max_tidx_ - 1`)
    virtual void UpdateSmoother(int tidx) = 0;  // Update smoother at time index
  };

  /// @brief Function signature for a custom fault-detection/outlier-rejection callback,
  /// registered via `KalmanFilter::SetFaultDetectionFunction`.
  typedef std::function<void(KalmanFilter*)> FaultDetectionFunction;

  /// @brief Extended Kalman Filter (EKF): linearized predict/update with Joseph-form
  /// covariance update, residual-based outlier rejection, and RTS smoothing.
  ///
  /// The standard recursive estimator used throughout LuPNT's navigation applications
  /// (e.g. LNSS receiver/satellite state estimation) -- registered with the `Filter` asset
  /// factory under the name `"EKF"` so it can be constructed from YAML config.
  class EKF : public KalmanFilter {
  protected:
    double outlier_threshold_ = 3.0;

    /// @brief Number of trailing state elements treated as Schmidt "consider"
    /// states (0 disables consider-state handling, i.e. a plain EKF).
    int n_consider_ = 0;

    std::vector<VecXd> x_prior_log_;  // Logged prior states
    std::vector<MatXd> P_prior_log_;  // Logged prior covariances
    std::vector<VecXd> x_pos_log_;    // Logged posterior states
    std::vector<MatXd> P_pos_log_;    // Logged posterior covariances
    std::vector<MatXd> stm_log_;      // Logged state transition matrices

  public:
    EKF() = default;

    /// @brief Construct an EKF from a YAML config node, reading the optional
    /// `outlier_threshold` field (default 3.0 sigma).
    EKF(Config& config);
    virtual ~EKF() = default;

    /// @brief Set the outlier-rejection threshold (in standard deviations of the
    /// normalized residual `dz_i / sqrt(S_ii)`) used by `RemoveOutliers`.
    ///
    /// @param outlier_threshold  Threshold in sigma; must be non-negative
    void SetOutlierThreshold(double outlier_threshold);

    /// @brief Configure the trailing `n_consider` elements of the state vector as
    /// Schmidt "consider" states, turning this `EKF` into a Schmidt (consider-
    /// parameter) Extended Kalman Filter (also available pre-configured as
    /// `SchmidtEKF`, see `schmidt_ekf.h`).
    ///
    /// Consider states are still propagated (`Predict`) and contribute to the
    /// measurement Jacobian `H_` and hence to the Kalman gain and covariance
    /// bookkeeping in `Update`, but are never corrected: `Update` zeros the
    /// trailing `n_consider` rows of the Kalman gain `K_` before applying the
    /// state correction, so their mean is left unchanged while their (co)variance
    /// -- including cross-covariance with the estimated states -- still updates
    /// consistently (the Joseph-form covariance update is valid for any gain,
    /// optimal or not). This is the standard "gain zeroing" implementation of the
    /// Schmidt-Kalman filter.
    ///
    /// @param n_consider  Number of trailing state elements to treat as consider
    ///                    states (0 disables consider-state handling)
    void SetConsiderStateCount(int n_consider);

    /// @brief Get the number of trailing consider states configured via
    /// `SetConsiderStateCount` (0 if this is a plain EKF).
    int GetConsiderStateCount() const { return n_consider_; }

    /// @brief Drop measurement components whose normalized residual
    /// `|dz_i| / sqrt(S_ii)` exceeds `outlier_threshold_`, shrinking `dz_`, `H_`, `R_`, and
    /// recomputing `S_` accordingly.
    ///
    /// Called by `Update` as the default fault-detection step (when no custom
    /// `f_fault_det_` is registered via `SetFaultDetectionFunction`), after the residual
    /// `dz_` and innovation covariance `S_` have been computed but before the Kalman gain
    /// is applied.
    void RemoveOutliers();

    /// @brief EKF predict step: propagate `x_`/`P_` via `f_dyn_` (linearized by the
    /// returned STM `F_`) and add the process noise `Q_` from `f_proc_`, i.e.
    /// `P_ = F_ P_ F_^T + Q_`.
    void Predict(Real t, const State* u = nullptr) override;

    /// @brief EKF update step: linearized measurement update with Joseph-form covariance
    /// propagation `P_ = (I - K H) P (I - K H)^T + K R K^T`, including outlier rejection
    /// via `RemoveOutliers` (or a registered fault-detection function) before the Kalman
    /// gain `K_ = P H^T S^-1` is applied.
    void Update(const VecX& z_true) override;

    /// @brief EKF implementation of `KalmanFilter::InitializeLogger`: allocates the
    /// `x_prior_log_`/`P_prior_log_`/`x_pos_log_`/`P_pos_log_`/`stm_log_`/`x_sm_`/`P_sm_`
    /// buffers for `max_tidx` epochs.
    void InitializeLogger(int max_tidx) override;

    /// @brief EKF implementation of `KalmanFilter::InitializeSmootherState`: seeds
    /// `x_sm_`/`P_sm_` at the final time index with the current `x_post_`/`P_post_`.
    void InitializeSmootherState() override;

    /// @brief EKF implementation of `KalmanFilter::LogFilterEstimate`: records the current
    /// prior/posterior state, covariance, and state-transition matrix `F_` at index `tidx`.
    void LogFilterEstimate(int tidx) override;

    /// @brief EKF implementation of `KalmanFilter::UpdateSmoother`: standard RTS smoother
    /// recursion using the logged prior/posterior states, covariances, and state-transition
    /// matrices.
    void UpdateSmoother(int tidx) override;
  };
}  // namespace lupnt
