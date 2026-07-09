/**
 * @file filter_utils.h
 * @author Stanford NAV Lab
 * @brief Utility functions for filters
 * @version 0.1
 * @date 2024-10-30
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include "lupnt/agents/agent.h"
#include "lupnt/devices/clock.h"
#include "lupnt/numerics/filters/filter.h"

namespace lupnt {

  /// @brief Wrap a Jacobian-free `MeasurementFunction` into a `FilterMeasurementFunction` by
  /// computing its Jacobian `H` via autodiff (`lupnt::jacobian`).
  ///
  /// Used by application setup code (e.g. `example_adaptive.cc`) to adapt a user-written
  /// measurement model `f_meas(x, R)` into the `FilterMeasurementFunction` signature
  /// expected by `Filter::SetMeasurementFunction`.
  ///
  /// @param f_meas  Jacobian-free measurement function `(x, R) -> z`
  /// @return        Equivalent `FilterMeasurementFunction` `(x, H, R) -> z` that also fills
  ///                in the measurement Jacobian `H` via autodiff
  FilterMeasurementFunction GetFilterMeasurementFunction(const MeasurementFunction& f_meas);

  /// @brief Wrap a Jacobian-free `DynamicsFunction` into a `FilterDynamicsFunction` by
  /// computing its state-transition Jacobian `F` via autodiff (`lupnt::jacobian`).
  ///
  /// Used by application setup code (e.g. `example_adaptive.cc`) to adapt a user-written
  /// dynamics model `f_dyn(x, t0, tf, u)` into the `FilterDynamicsFunction` signature
  /// expected by `Filter::SetDynamicsFunction`.
  ///
  /// @param f_dyn  Jacobian-free dynamics function `(x, t0, tf, u) -> x_tf`
  /// @return       Equivalent `FilterDynamicsFunction` `(x, t0, tf, u, F) -> x_tf` that also
  ///               fills in the state-transition Jacobian `F` via autodiff
  FilterDynamicsFunction GetFilterDynamicsFunction(const DynamicsFunction& f_dyn);

  /// @brief Identity pass-through for a `ProcessNoiseFunction` (returns `f_proc` unchanged).
  ProcessNoiseFunction GetProcessNoiseFunction(const ProcessNoiseFunction& f_proc);

  /// @brief Build a block-diagonal initial covariance for a position/velocity/clock-bias/
  /// clock-drift state: `diag(sigma_r^2 * I3, sigma_v^2 * I3, sigma_b^2, sigma_d^2)`.
  ///
  /// Used by application setup code to initialize a filter's `P0` (via
  /// `Filter::SetCovariance`) for the common 8-element `[r(3), v(3), clock_bias,
  /// clock_drift]` state.
  ///
  /// @param sigma_r  Initial position 1-sigma uncertainty [m]
  /// @param sigma_v  Initial velocity 1-sigma uncertainty [m/s]
  /// @param sigma_b  Initial clock-bias 1-sigma uncertainty [s]
  /// @param sigma_d  Initial clock-drift 1-sigma uncertainty [s/s]
  /// @return         Initial covariance matrix, size [8 x 8]
  MatXd InitialCovariancePosVelClock(double sigma_r, double sigma_v, double sigma_b,
                                     double sigma_d);

  /// @brief Compute the position/velocity process-noise coefficients `[C_11, C_21, C_22]`
  /// for a piecewise-constant white-acceleration (state-noise-compensation) model over a
  /// time step `dt`.
  ///
  /// Used by `ProcessNoisePosVel` to build the position/velocity process noise block, and
  /// directly by `ProcessNoise::ComputeAccNoise` (ASNC algorithm) and
  /// `ProcessNoise::SolveForQtilde` to relate the acceleration noise PSD to the
  /// observed position/velocity covariance growth.
  ///
  /// @param dt  Time step [s] (must be positive)
  /// @return    Coefficients `[C_11, C_21, C_22]` for the position/velocity process-noise
  ///            blocks
  Vec3d ProcessNoisePosVelCoeffs(Real dt);

  /// @brief Compute the position/velocity/acceleration process-noise coefficients
  /// `[C_11, C_21, C_22, C_31, C_32, C_33]` (one column per acceleration component) for a
  /// first-order Gauss-Markov (exponentially-correlated) acceleration model with time
  /// constants `1/beta` over a time step `dt`.
  ///
  /// Used by `ProcessNoisePosVelAcc` to build the position/velocity/acceleration process
  /// noise block, and directly by `ProcessNoise::ComputeAccNoise` (ADMC algorithm) and
  /// `ProcessNoise::SolveForQtilde`.
  ///
  /// @param dt    Time step [s] (must be positive)
  /// @param beta  Gauss-Markov inverse time constants per acceleration component, size [n]
  /// @return      Coefficient matrix, size [6 x n] (rows: C_11, C_21, C_22, C_31, C_32, C_33)
  Mat6Xd ProcessNoisePosVelAccCoeffs(Real dt, const VecXd& beta);

  /// @brief Build the position/velocity process-noise covariance block for a
  /// piecewise-constant white-acceleration model, from the per-axis acceleration noise
  /// `Q_a` and `ProcessNoisePosVelCoeffs(dt)`.
  ///
  /// Called by `ProcessNoise::ComputeProcessNoise` (SNC/ASNC algorithms) and
  /// `ProcessNoisePosVelClock` to form the `[r, v]` block of the EKF/UKF process-noise
  /// matrix `Q` passed to `Filter::SetProcessNoiseFunction`.
  ///
  /// @param Q_a  Per-axis acceleration noise: either a diagonal vector or a full matrix,
  ///             size [n] or [n x n]
  /// @param dt   Time step [s]
  /// @return     Position/velocity process noise covariance, size [2n x 2n]
  MatXd ProcessNoisePosVel(const MatXd& Q_a, Real dt);

  /// @brief Build the position/velocity/acceleration process-noise covariance block for a
  /// first-order Gauss-Markov acceleration model, from the per-axis driving noise `Q_a`,
  /// time step `dt`, and inverse time constants `beta`.
  ///
  /// Called by `ProcessNoise::ComputeProcessNoise` (ADMC algorithm) to form the
  /// `[r, v, a]` block of the process-noise matrix `Q`.
  ///
  /// @param Q_a   Per-axis driving-noise: either a diagonal vector or a full matrix,
  ///              size [n] or [n x n]
  /// @param dt    Time step [s]
  /// @param beta  Gauss-Markov inverse time constants per acceleration component, size [n]
  /// @return      Position/velocity/acceleration process noise covariance, size [3n x 3n]
  MatXd ProcessNoisePosVelAcc(const MatXd& Q_a, Real dt, const VecXd& beta);

  /// @brief Build the clock-state process-noise covariance block for a given clock model
  /// and time step, via `ClockDynamics::TwoStateNoise`/`ThreeStateNoise`.
  ///
  /// Called by `ProcessNoisePosVelClock` to form the clock block of the combined
  /// position/velocity/clock process-noise matrix `Q`.
  ///
  /// @param clock_model  Clock noise model (e.g. OCXO, USO, CSAC, RAFS)
  /// @param n_clk        Clock state dimension: 2 (bias, drift) or 3 (bias, drift,
  ///                      drift-rate)
  /// @param dt           Time step [s]
  /// @return             Clock process noise covariance, size [n_clk x n_clk]
  MatXd ProcessNoiseClock(ClockModel cmodel, int n_clk, Real dt);

  /// @brief Build a block-diagonal process-noise covariance for `n_sat` independent
  /// position/velocity/clock-bias/clock-drift states, combining `ProcessNoisePosVel` (with
  /// an isotropic per-axis acceleration noise `sigma_a`) and `ProcessNoiseClock` per
  /// satellite.
  ///
  /// Intended as the `Q` returned by a `ProcessNoiseFunction` registered via
  /// `Filter::SetProcessNoiseFunction` for a multi-satellite position/velocity/clock filter
  /// state.
  ///
  /// @param cmodel  Clock noise model for each satellite's clock
  /// @param n_clk   Per-satellite block size (position/velocity + clock dimensions, e.g. 8
  ///                for `[r(3), v(3), bias, drift]`)
  /// @param sigma_a Per-axis acceleration noise 1-sigma [m/s^2]
  /// @param n_sat   Number of satellites/states (default 1)
  /// @return        Block-diagonal process noise covariance, size [n_clk*n_sat x n_clk*n_sat]
  MatXd ProcessNoisePosVelClock(ClockModel cmodel, int n_clk, double sigma_a, int n_sat = 1);

  /// @brief Build the constant-velocity position/velocity state-transition matrix
  /// `Phi = [[I, dt*I], [0, I]]` for `n`-dimensional position.
  ///
  /// Used as the dynamics STM `F` (e.g. returned from a `FilterDynamicsFunction`) for
  /// simple constant-velocity filter models, such as in `example_adaptive.cc`.
  ///
  /// @param dt  Time step [s]
  /// @param n   Dimension of the position (and velocity) sub-state
  /// @return    State-transition matrix, size [2n x 2n]
  MatXd StateTransitionMatrixPosVel(double dt, int n);

  /// @brief Build the position/velocity/acceleration state-transition matrix for a
  /// first-order Gauss-Markov acceleration model with inverse time constants `beta`,
  /// combining the constant-velocity block (`StateTransitionMatrixPosVel`) with the
  /// exponential acceleration decay and its position/velocity coupling terms.
  ///
  /// Used as the dynamics STM `F` for filter models that include a Gauss-Markov
  /// acceleration state (ADMC), such as in `example_adaptive.cc`.
  ///
  /// @param dt    Time step [s]
  /// @param beta  Gauss-Markov inverse time constants per acceleration component, size [n]
  /// @return      State-transition matrix, size [3n x 3n]
  MatXd StateTransitionMatrixPosVelAcc(double dt, const VecXd& beta);

  /// @brief Concatenate the 8-element `[r(3), v(3), clock_bias, clock_drift]` true state
  /// vectors of a set of satellites into a single stacked vector.
  ///
  /// Used to build the `true_state` reference against which a multi-satellite filter's
  /// estimate is compared (e.g. via `ComputeEstimationErrorPVC`).
  ///
  /// @param sats  Satellites whose `GetState()` vectors are concatenated
  /// @return      Stacked true-state vector, size [8 * sats.size()]
  VecX ConstructTrueStateVecFromSats(const std::vector<Ptr<Agent>>& sats);

  /// @brief Compute the position/velocity/clock-bias/clock-drift estimation errors for one
  /// satellite, comparing its true state (`sat->GetState()`) against the filter's posterior
  /// estimate starting at index `start_idx`.
  ///
  /// Used by application logging/printing code (e.g. with `PrintEKFProgressPVC`) to report
  /// per-satellite navigation accuracy at each filter step.
  ///
  /// @param sat        Satellite providing the true state
  /// @param filter     Filter providing the posterior state estimate (`GetStatePost`)
  /// @param start_idx  Index of this satellite's position sub-state within the filter's
  ///                    state vector
  /// @return           Errors `[pos_err [mm], vel_err [um/s], clk_bias_err [m],
  ///                    clk_drift_err [m/s]]` (scaled by the constants used in the
  ///                    implementation)
  VecXd ComputeEstimationErrorPVC(const Ptr<Agent>& sat, Filter* filter, int start_idx);

  /// @brief Overload of `ComputeEstimationErrorPVC` for multiple satellites: stacks the
  /// per-satellite 4-element error vectors, assuming each satellite occupies an 8-element
  /// block of the filter state.
  ///
  /// @param sats    Satellites whose true states are compared against the filter estimate
  /// @param filter  Filter providing the posterior state estimate
  /// @return        Stacked error vector, size [4 * sats.size()]
  VecXd ComputeEstimationErrorPVC(const std::vector<Ptr<Agent>>& sats, Filter* filter);

}  // namespace lupnt
