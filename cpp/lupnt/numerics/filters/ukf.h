/**
 * @file ukf.h
 * @author Stanford NAV LAB
 * @brief  Unscented Kalman Filter
 * @version 0.1
 * @date 2024-12-30
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include "lupnt/core/definitions.h"
#include "lupnt/numerics/filters/filter.h"

namespace lupnt {

  /// @brief Unscented Kalman Filter (UKF): predict/update via the unscented transform
  /// (sigma-point propagation) instead of linearization.
  ///
  /// An alternative to `EKF` for nonlinear dynamics/measurement models where Jacobians are
  /// unavailable or linearization error is too large; constructed and driven the same way
  /// as other `Filter` subclasses (`SetDynamicsFunction`/`SetProcessNoiseFunction`/
  /// `SetMeasurementFunction`, then `Predict`/`Update`).
  class UKF : public Filter {
  private:
    // Unscented transform parameters
    double alpha_ = 1e-3;
    double beta_ = 2.0;
    double kappa_ = 0.0;

    // Derived quantities
    int n_x_ = 0;        // dimension of the state
    int n_sigma_ = 0;    // number of sigma points = 2*n_x_ + 1
    double lambda_ = 0;  // alpha_^2*(n_x_+kappa_) - n_x_

    // Weights for mean and covariance
    VecXd w_m_;  // mean weights
    VecXd w_c_;  // covariance weights

    /// @brief Recompute the unscented-transform scaling factor `lambda_`, sigma-point count
    /// `n_sigma_ = 2*n_x_ + 1`, and mean/covariance sigma-point weights `w_m_`/`w_c_` from
    /// `alpha_`, `beta_`, `kappa_`, and `n_x_`.
    ///
    /// Called by `Initialize` (once `n_x_` is known) and by `SetAlpha`/`SetBeta`/`SetKappa`
    /// whenever a unscented-transform tuning parameter changes, and lazily by `Predict` on
    /// first use.
    void InitializeUkfParams();

  protected:
    MatXd F_;  // State transition matrix
    MatXd H_;  // Measurement matrix
    MatXd Q_;  // Process noise cov
    MatXd R_;  // Measurement noise cov

    VecX dy_;       // Measurement residual
    VecX dx_;       // State update
    VecX z_true_;   // Observed measurement
    VecX z_prior_;  // Predicted measurement

    MatXd S_;  // Innovation cov
    MatXd K_;  // Kalman gain

  public:
    UKF() = default;

    /// @brief Initialize the UKF's epoch, state, and covariance (including prior/posterior
    /// copies), set the state dimension `n_x_` from `x0`, and (re)compute the unscented
    /// transform weights via `InitializeUkfParams`.
    ///
    /// Called once by application setup code before the first `Predict`/`Update` call.
    ///
    /// @param t0  Initial epoch [s]
    /// @param x0  Initial state estimate, size [n_x]
    /// @param P0  Initial state covariance, size [n_x x n_x]
    void Initialize(const Real t0, const State& x0, const MatXd& P0) {
      t_ = t0;
      x_ = x0;
      P_ = P0;

      n_x_ = x0.size();
      x_prior_ = x0;
      P_prior_ = P0;
      x_post_ = x0;
      P_post_ = P0;

      InitializeUkfParams();
    };

    /// @brief Generate the `2*n_x_ + 1` unscented-transform sigma points for a given mean
    /// `state` and covariance `cov`: the central point `state`, plus `state +/-
    /// chol((n_x_ + lambda_) * cov)` columns.
    ///
    /// Called by `Predict` (on the prior `x_`/`P_`, propagated through `f_dyn_`) and by
    /// `Update` (on the predicted `x_`/`P_`, passed through `f_meas_`).
    ///
    /// @param state  Mean state about which to spread the sigma points, size [n_x]
    /// @param cov    Covariance matrix, size [n_x x n_x]
    /// @return       Sigma points, size [n_x x n_sigma_] (n_sigma_ = 2*n_x_ + 1)
    MatX ComputeSigmaPoints(const State& state, const MatXd& cov);

    /// @brief Recombine a set of (already-propagated) sigma points into a mean and
    /// covariance using the unscented-transform weights `w_m_`/`w_c_`, optionally adding a
    /// process noise covariance `Q`.
    ///
    /// Called by `Predict` to recombine the propagated dynamics sigma points (with `Q` =
    /// the process noise from `f_proc_`) into the predicted state mean/covariance.
    ///
    /// @param sigma_points  Propagated sigma points, size [n_x_ x n_sigma_]
    /// @param cov_out       Output covariance of the recombined state, size [n_x_ x n_x_]
    /// @param Q             Optional process noise covariance to add, size [n_x_ x n_x_]
    /// @return              Mean of the sigma points, size [n_x_]
    State UnscentedTransform(const MatX& sigma_points, MatXd& cov_out, const MatXd* Q = nullptr);

    /// @brief Set the unscented-transform spread parameter `alpha` (typically small,
    /// e.g. 1e-3) and recompute the sigma-point weights.
    void SetAlpha(double alpha) {
      alpha_ = alpha;
      InitializeUkfParams();
    }
    /// @brief Set the unscented-transform secondary scaling parameter `beta` (2 is optimal
    /// for Gaussian distributions) and recompute the sigma-point weights.
    void SetBeta(double beta) {
      beta_ = beta;
      InitializeUkfParams();
    }
    /// @brief Set the unscented-transform secondary scaling parameter `kappa` (typically 0
    /// or `3 - n_x`) and recompute the sigma-point weights.
    void SetKappa(double kappa) {
      kappa_ = kappa;
      InitializeUkfParams();
    }

    /// @brief UKF predict step: generate sigma points from `(x_, P_)` via
    /// `ComputeSigmaPoints`, propagate each through `f_dyn_` (no Jacobian needed), then
    /// recombine via `UnscentedTransform` (adding process noise `Q_` from `f_proc_`) to
    /// form the predicted `x_prior_`/`P_prior_`, which also become the new `x_`/`P_`.
    void Predict(Real t, const State* u = nullptr) override;

    /// @brief UKF update step: generate sigma points from the predicted `(x_, P_)`,
    /// transform each through `f_meas_` to form the measurement sigma points, recombine
    /// them (with weights `w_m_`/`w_c_`) into the predicted measurement mean and innovation
    /// covariance `S_`, compute the state-measurement cross-covariance, form the Kalman
    /// gain `K_ = cross_cov * S_^-1`, and apply the correction `dx_ = K_ * dy_` to update
    /// `x_post_`/`P_post_` (which also become the new `x_`/`P_`).
    void Update(const VecX& z_true) override;

    // Interface
    /// @brief Get the measurement residual `dy = z_true - meas_mean` from the last `Update`
    /// call, size [n_z].
    VecXd GetMeasurementResidual() { return dy_.cast<double>(); }
    /// @brief Get the Kalman gain `K_ = cross_cov * S_^-1` from the last `Update` call,
    /// size [n_x x n_z].
    MatXd GetKalmanGain() { return K_; }
    /// @brief Get the measurement noise covariance `R_` (from the first sigma point's
    /// `f_meas_` evaluation) from the last `Update` call, size [n_z x n_z].
    MatXd GetMeasurementNoiseCov() { return R_; }
    /// @brief Get `H_` (unused by the UKF's sigma-point update; retained for interface
    /// compatibility with `KalmanFilter`).
    MatXd GetMeasurementJacobian() { return H_; }
    /// @brief Get the current measurement dimension (number of rows of `H_`; see
    /// `GetMeasurementJacobian`).
    int GetMeasurementSize() { return H_.rows(); }
    /// @brief Get `Q_` (unused by the UKF's sigma-point predict step, which adds process
    /// noise directly in `UnscentedTransform`; retained for interface compatibility with
    /// `KalmanFilter`).
    MatXd GetProcessNoise() { return Q_; }
    /// @brief Get `F_` (unused by the UKF's sigma-point predict step, which has no explicit
    /// state-transition Jacobian; retained for interface compatibility with
    /// `KalmanFilter`).
    MatXd GetStateJacobian() { return F_; }
    /// @brief Get the innovation covariance `S_` from the last `Update` call,
    /// size [n_z x n_z].
    MatXd GetInnovationCov() { return S_; }
    /// @brief Same as `GetMeasurementNoiseCov`; get the measurement noise covariance `R_`.
    MatXd GetMeasurementCov() { return R_; }
    /// @brief Get the state correction `dx_ = K_ * dy_` from the last `Update` call,
    /// size [n_x].
    VecXd GetStateCorrection() { return dx_.cast<double>(); }
    /// @brief Get the observed measurement vector from the last `Update` call, size [n_z].
    VecXd GetTrueMeasurement() { return z_true_.cast<double>(); }
    /// @brief Get the predicted measurement mean (recombined from the measurement sigma
    /// points) from the last `Update` call, size [n_z].
    VecXd GetPredictedMeasurement() { return z_prior_.cast<double>(); }
  };

}  // namespace lupnt
