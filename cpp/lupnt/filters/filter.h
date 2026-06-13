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
#include <functional>

#include "lupnt/core/config.h"
#include "lupnt/core/data_logger.h"
#include "lupnt/core/object.h"
#include "lupnt/states/state.h"

namespace lupnt {

  // Dynamics and Measurement Function

  /// @brief Function signature for propagating the filter state from `t0` to `tf` and
  /// returning the state-transition (Jacobian) matrix `F` of the propagation.
  ///
  /// Stored in `Filter::f_dyn_` (via `SetDynamicsFunction`) and called once per `Predict`
  /// step by `EKF::Predict`, `UDUEKF::Predict`, and `UKF::Predict` (the latter with
  /// `F = nullptr` since the UKF propagates sigma points instead of a Jacobian). Typically
  /// produced from a plain `DynamicsFunction` via `GetFilterDynamicsFunction`.
  ///
  /// @param x  State at time `t0`
  /// @param t0 Current time [s]
  /// @param tf End time [s]
  /// @param u  Optional control/input state (nullptr if unused)
  /// @param F  Output state-transition (Jacobian) matrix d(x_tf)/d(x_t0), size [n_x x n_x]
  /// @return   Propagated state at time `tf`
  typedef std::function<State(const State& x, Real t0, Real tf, const State* u, MatXd* F)>
      FilterDynamicsFunction;

  /// @brief Function signature for propagating the filter state from `t0` to `tf` without
  /// computing a Jacobian.
  ///
  /// Wrapped into a `FilterDynamicsFunction` (which additionally produces the Jacobian via
  /// autodiff) by `GetFilterDynamicsFunction`.
  ///
  /// @param x  State at time `t0`
  /// @param t0 Current time [s]
  /// @param tf End time [s]
  /// @param u  Optional control/input state (nullptr if unused)
  /// @return   Propagated state at time `tf`
  typedef std::function<State(const State& x, Real t0, Real tf, const State* u)> DynamicsFunction;

  /// @brief Function signature for computing the process noise covariance accumulated over
  /// `[t0, tf]`.
  ///
  /// Stored in `Filter::f_proc_` (via `SetProcessNoiseFunction`) and called once per
  /// `Predict` step by `EKF::Predict`, `UDUEKF::Predict`, and `UKF::Predict` to form `Q_`,
  /// which is added to the propagated state covariance.
  ///
  /// @param x  State at time `t0`
  /// @param t0 Current time [s]
  /// @param tf End time [s]
  /// @return   Process noise covariance matrix, size [n_x x n_x]
  typedef std::function<MatXd(const State& x, Real t0, Real tf)> ProcessNoiseFunction;

  /// @brief Function signature for predicting a measurement from the current state, and
  /// computing its Jacobian `H` and noise covariance `R`.
  ///
  /// Stored in `Filter::f_meas_` (via `SetMeasurementFunction`) and called once per `Update`
  /// step by `EKF::Update`/`UDUEKF::Update` (and once per sigma point plus once for sizing by
  /// `UKF::Update`) to form the measurement residual `dz = z_true - z_prior`. Typically
  /// produced from a plain `MeasurementFunction` via `GetFilterMeasurementFunction`, and
  /// often built from a measurement model's `CreateFunction` (e.g.
  /// `GnssMeasurement::CreateFunction`).
  ///
  /// @param x  Current state estimate
  /// @param H  Output measurement Jacobian d(z)/d(x), size [n_z x n_x]
  /// @param R  Output measurement noise covariance, size [n_z x n_z]
  /// @return   Predicted measurement vector, size [n_z]
  typedef std::function<VecXd(const State& x, MatXd* H, MatXd* R)> FilterMeasurementFunction;

  /// @brief Function signature for predicting a measurement from the current state and
  /// computing its noise covariance `R`, without a Jacobian.
  ///
  /// Wrapped into a `FilterMeasurementFunction` (which additionally produces the Jacobian
  /// via autodiff) by `GetFilterMeasurementFunction`.
  ///
  /// @param x  Current state estimate
  /// @param R  Output measurement noise covariance, size [n_z x n_z]
  /// @return   Predicted measurement vector, size [n_z]
  typedef std::function<VecXd(const State& x, MatXd* R)> MeasurementFunction;

  /// @brief Abstract base class for recursive (sequential) navigation filters that maintain
  /// a state estimate `x_` and covariance `P_`, propagated by `f_dyn_`/`f_proc_` and
  /// corrected by `f_meas_`.
  ///
  /// Concrete filters (EKF, UDUEKF, UKF) implement the `Predict`/`Update` pair below with
  /// algorithm-specific propagation/correction math. Application code (e.g. the LNSS/rover
  /// applications) configures a `Filter` via `SetDynamicsFunction`/`SetProcessNoiseFunction`/
  /// `SetMeasurementFunction`, initializes `x_`/`P_` via `SetState`/`SetCovariance`, and then
  /// drives the estimation loop by alternating `Predict(t)` (propagate to time `t`) and
  /// `Update(z_true)` (incorporate a new measurement vector).
  class Filter : public Object<Filter>, public DataLogger {
  protected:
    Config config_;
    std::string name_;

    FilterDynamicsFunction f_dyn_;
    ProcessNoiseFunction f_proc_;
    FilterMeasurementFunction f_meas_;

    Real t_ = 0.0;
    State x_;
    State x_prior_;
    State x_post_;

    MatXd P_;
    MatXd P_prior_;
    MatXd P_post_;
    MatXd P_bar_;

  public:
    Filter() = default;

    /// @brief Construct a filter from a YAML config node, taking its `name` field (or the
    /// auto-generated object id if absent).
    Filter(Config& config) : config_(config) {
      name_ = config["name"] ? config["name"].as<std::string>() : GetId();
    }
    virtual ~Filter() = default;

    /// @brief Set the filter's display/log name (used as the prefix for `Log` entries).
    void SetName(const std::string& name) { name_ = name; }
    /// @brief Get the filter's display/log name.
    std::string GetName() const { return name_; }

    /// @brief Register the state-propagation function used by `Predict`.
    void SetDynamicsFunction(FilterDynamicsFunction f_dyn) { f_dyn_ = f_dyn; }
    /// @brief Register the process-noise covariance function used by `Predict`.
    void SetProcessNoiseFunction(ProcessNoiseFunction f_proc) { f_proc_ = f_proc; }
    /// @brief Register the measurement-prediction function used by `Update`.
    void SetMeasurementFunction(FilterMeasurementFunction f_meas) { f_meas_ = f_meas; }

    /// @brief Set the filter's current epoch `t_` (used as the `t0` argument of the next
    /// `Predict` call).
    void SetTime(Real t) { t_ = t; }

    /// @brief Initialize the state estimate `x_` (and prior/posterior copies) at filter setup.
    ///
    /// Called once by application setup code before the estimation loop begins. Overridden
    /// by `UDUEKF::SetCovariance`'s sibling (see `SetCovariance`) but `SetState` itself is
    /// not specialized by the UD filter.
    ///
    /// @param x  Initial state estimate
    virtual void SetState(const State& x) {
      x_ = x;
      x_prior_ = x;
      x_post_ = x;
    }

    /// @brief Initialize the state covariance `P_` (and prior/posterior copies) at filter
    /// setup.
    ///
    /// Called once by application setup code before the estimation loop begins.
    /// `UDUEKF::SetCovariance` overrides this to additionally factor `P` into `U_`/`D_diag_`.
    ///
    /// @param P  Initial state covariance matrix, size [n_x x n_x]
    virtual void SetCovariance(const MatXd& P) {
      P_ = P;
      P_prior_ = P;
      P_post_ = P;
    }

    /// @brief Propagate the state estimate and covariance from the current epoch `t_` to
    /// time `t`, using `f_dyn_`/`f_proc_`, storing the result in `x_`/`P_` (and the prior
    /// copies `x_prior_`/`P_prior_`).
    ///
    /// Called once per estimation step (before `Update`) by the application's filtering
    /// loop, e.g. once per epoch in the LNSS/rover navigation applications. `EKF::Predict`
    /// linearizes via the STM `F_` from `f_dyn_`; `UKF::Predict` propagates unscented
    /// transform sigma points; `UDUEKF::Predict` does the same as `EKF::Predict` but
    /// maintains `P_` as a UDU (Bierman-Thornton) factorization.
    ///
    /// @param t  Target epoch to propagate to [s]
    /// @param u  Optional control/input state passed through to `f_dyn_` (nullptr if unused)
    virtual void Predict(Real t, const State* u = nullptr) = 0;

    /// @brief Incorporate a new measurement vector `z_true` into the state estimate and
    /// covariance via `f_meas_`, updating `x_`/`P_` (and the posterior copies
    /// `x_post_`/`P_post_`).
    ///
    /// Called once per estimation step (after `Predict`) whenever a new measurement is
    /// available. `EKF::Update` performs a linearized Joseph-form Kalman update;
    /// `UKF::Update` performs an unscented-transform-based update; `UDUEKF::Update` performs
    /// a sequential Carlson (UD/Bierman-Thornton) update.
    ///
    /// @param z_true  Observed measurement vector, size [n_z]
    virtual void Update(const VecX& z_true) = 0;

    /// @brief Get the filter's current epoch [s].
    Real GetTime() { return t_; }
    /// @brief Get the current state estimate.
    State GetState() { return x_; }
    /// @brief Get the predicted (pre-update) state estimate from the last `Predict` call.
    State GetStatePrior() { return x_prior_; }
    /// @brief Get the corrected (post-update) state estimate from the last `Update` call.
    State GetStatePost() { return x_post_; }

    /// @brief Get the current state covariance, size [n_x x n_x].
    MatXd GetCovariance() { return P_; }
    /// @brief Get the predicted (pre-update) state covariance from the last `Predict` call.
    MatXd GetCovariancePrior() { return P_prior_; }
    /// @brief Get the corrected (post-update) state covariance from the last `Update` call.
    MatXd GetCovariancePost() { return P_post_; }
    /// @brief Get the propagated-only covariance `F P F^T` (before adding process noise)
    /// from the last `Predict` call.
    MatXd GetCovarianceBar() { return P_bar_; }

    /// @brief Compute the (unweighted-by-size) chi-square measurement residual sum
    /// `sum_i dz_i^2 / R_ii` for `z_true` against the measurement predicted at `x_est`,
    /// skipping any NaN or non-positive-variance components.
    ///
    /// Used as a goodness-of-fit / outlier diagnostic; returns 0 if `z_true` is empty.
    ///
    /// @param z_true  Observed measurement vector, size [n_z] (or empty)
    /// @param x_est   State at which to evaluate the measurement function `f_meas_`
    /// @return        Sum of squared, variance-normalized measurement residuals
    double ComputeResidualRMS(VecXd z_true, const State& x_est);

    /// @brief Log the current epoch, state, and covariance diagonals (current,
    /// prior, and posterior) to the `DataLogger` under the `<name_>/...` keys.
    ///
    /// Called once per estimation step by the application's logging loop to record filter
    /// history for post-run analysis/plotting.
    ///
    /// @param time  Current simulation time [s] (unused; `t_` is logged instead)
    virtual void Log(Real time) override;
  };
}  // namespace lupnt
