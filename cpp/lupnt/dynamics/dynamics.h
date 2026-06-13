/**
 * @file dynamics.h
 * @author Stanford NAV LAB
 * @brief Interface for Dynamics
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <yaml-cpp/yaml.h>

#include <map>

#include "lupnt/core/config.h"
#include "lupnt/states/params.h"
#include "lupnt/states/state.h"

namespace lupnt {

  /// @brief Base interface for all state-propagation models in LuPNT (orbit, attitude,
  /// clock, IMU, surface, parameter, and joint dynamics).
  ///
  /// A `Dynamics` subclass owns the "physics" for one state vector: it defines how to
  /// advance a `State` from `t0` to `tf` (`Propagate`), optionally returning a state
  /// transition matrix (STM) via automatic differentiation, and exposes a
  /// `ParamState` of force-model/clock/etc. parameters that filters
  /// (filters/) can read, perturb, or estimate. Numerical subclasses
  /// (NumericalDynamics and its descendants) additionally provide a `ComputeRates`
  /// right-hand-side consumed by the integrators in numerics/integrator.h.
  class Dynamics {
  protected:
    bool print_progress_ = false;  ///< Flag to indicate whether to print progress.
    ParamState params_;            ///< Parameters for the dynamics, initialized as empty.

  public:
    Dynamics() = default;

    /// @brief Construct from a YAML configuration node.
    ///
    /// Reads the optional `print_progress` flag. Called by the
    /// `AssetFactory<Dynamics, Config&>`-based constructors of concrete dynamics
    /// classes (e.g. NumericalDynamics, ClockDynamics, AttitudeDynamics) when an
    /// agent/asset is built from a simulation YAML config.
    ///
    /// @param config  YAML configuration node for this dynamics model.
    Dynamics(Config& config);
    virtual ~Dynamics() = default;

    /// @brief Enable/disable printing a progress bar during multi-step Propagate calls.
    void SetPrintProgress(bool print) { print_progress_ = print; }
    /// @brief Return whether progress printing is enabled.
    bool GetPrintProgress() { return print_progress_; }

    /// @brief Propagate the state from `t0` to `tf` without computing a state
    /// transition matrix.
    ///
    /// This is the core propagation entry point implemented by every concrete
    /// dynamics class (e.g. NumericalDynamics::Propagate,
    /// AnalyticalDynamics::Propagate, ClockDynamics::Propagate). It is called once
    /// per integration/measurement step by simulation drivers in applications/ and
    /// by agents/ to advance an object's true or estimated state.
    ///
    /// @param x0  Initial state at time `t0` (layout/units defined by the subclass).
    /// @param t0  Initial epoch [s, TDB since J2000 for orbit/attitude dynamics].
    /// @param tf  Final epoch [s, TDB since J2000 for orbit/attitude dynamics].
    /// @param u   Optional control/forcing input (e.g. thrust, RTN reference state);
    ///            nullptr if unused.
    /// @return    Propagated state at time `tf`.
    virtual State Propagate(const State& x0, Real t0, Real tf, const State* u = nullptr) = 0;

    /// @brief Propagate the state and compute the state transition matrix (STM)
    /// d(xf)/d(x0) via automatic differentiation.
    ///
    /// Default implementation: if `stm == nullptr`, forwards to the no-STM
    /// `Propagate` overload; otherwise wraps `Propagate(x0, t0, tf, u)` in an
    /// autodiff Jacobian (`jacobian(...)`) with respect to `x0`. Used by filters
    /// (filters/) that need the linearized state transition for covariance
    /// propagation (e.g. EKF/UKF time updates).
    ///
    /// @param x0   Initial state at time `t0`.
    /// @param t0   Initial epoch [s, TDB since J2000].
    /// @param tf   Final epoch [s, TDB since J2000].
    /// @param u    Optional control/forcing input; nullptr if unused.
    /// @param stm  Output: state transition matrix d(xf)/d(x0); left unmodified if
    ///             nullptr.
    /// @return     Propagated state at time `tf`.
    virtual State Propagate(const State& x0, Real t0, Real tf, const State* u, MatXd* stm);

    /// @brief Propagate the state to a sequence of output epochs.
    ///
    /// Repeatedly calls `Propagate(x0, t0, tf, u)` between consecutive entries of
    /// `tfs`, accumulating each result as a row of the returned matrix. Used by
    /// applications/ and analysis scripts to generate a full trajectory time
    /// history in one call, with optional progress-bar output (see
    /// SetPrintProgress).
    ///
    /// @param x0   Initial state at `tfs(0)`.
    /// @param tfs  Vector of output epochs [s, TDB since J2000], including the
    ///             initial epoch as `tfs(0)`.
    /// @param u    Optional control/forcing input applied at every step; nullptr if
    ///             unused.
    /// @return     Matrix whose i-th row is the propagated state at `tfs(i)`.
    virtual MatX Propagate(const State& x0, const VecX& tfs, const State* u = nullptr);

    /// @brief Propagate the state after first installing a given parameter set.
    ///
    /// Calls `SetParams(params)` then `Propagate(x0, t0, tf, u)`. Used by filters
    /// and consider-covariance analyses that need to evaluate the dynamics under a
    /// specific (e.g. perturbed or estimated) parameter vector without permanently
    /// mutating the dynamics object's stored parameters beforehand.
    ///
    /// @param x0      Initial state at time `t0`.
    /// @param t0      Initial epoch [s, TDB since J2000].
    /// @param tf      Final epoch [s, TDB since J2000].
    /// @param params  Parameter values/names to install before propagating (see
    ///                 ParamState, GetParams/SetParams).
    /// @param u       Optional control/forcing input; nullptr if unused.
    /// @return        Propagated state at time `tf`.
    virtual State PropagateWithParams(const State& x0, Real t0, Real tf, const ParamState& params,
                                      const State* u = nullptr);

    /// @brief Propagate the state and parameters, returning both the state-transition
    /// and parameter-sensitivity Jacobians.
    ///
    /// If both `stm_state` and `stm_param` are nullptr, forwards to the no-Jacobian
    /// `PropagateWithParams` overload. Otherwise computes d(xf)/d(x0) into
    /// `stm_state` and d(xf)/d(params) into `stm_param` via autodiff. Used by
    /// batch least-squares and EKF/UKF filters (filters/) that jointly estimate the
    /// state and dynamics parameters (e.g. SRP/drag coefficients, clock-bias
    /// parameters).
    ///
    /// @param x0        Initial state at time `t0`.
    /// @param t0        Initial epoch [s, TDB since J2000].
    /// @param tf        Final epoch [s, TDB since J2000].
    /// @param params    Parameter values/names to install before propagating.
    /// @param u         Optional control/forcing input; nullptr if unused.
    /// @param stm_state State-transition Jacobian d(xf)/d(x0); unmodified if nullptr.
    /// @param stm_param Parameter-sensitivity Jacobian d(xf)/d(params); unmodified if
    ///                   nullptr.
    /// @return          Propagated state at time `tf`.
    virtual State PropagateWithParams(const State& x0, Real t0, Real tf, const ParamState& params,
                                      const State* u, MatXd* stm_state, MatXd* stm_param);

    /// @brief Propagate the state with the given parameters to a sequence of output
    /// epochs.
    ///
    /// Sequence-output counterpart of PropagateWithParams, analogous to
    /// `Propagate(x0, tfs, u)` but holding `params` fixed across all steps.
    ///
    /// @param x0      Initial state at `tfs(0)`.
    /// @param tfs     Vector of output epochs [s, TDB since J2000].
    /// @param params  Parameter values/names to install before propagating.
    /// @param u       Optional control/forcing input applied at every step; nullptr
    ///                 if unused.
    /// @return        Matrix whose i-th row is the propagated state at `tfs(i)`.
    virtual MatX PropagateWithParams(const State& x0, const VecX& tfs, const ParamState& params,
                                     const State* u = nullptr);

    /// @brief Replace the dynamics' parameter vector (e.g. SRP/drag coefficients,
    /// clock model parameters) wholesale.
    void SetParams(const ParamState& params) { params_ = params; }

    /// @brief Set a single named parameter's value in-place.
    ///
    /// Looks up `key` among the current parameter names and overwrites that entry's
    /// value with `value`. Used by force-model setters (e.g.
    /// NBodyDynamics::SetSrpCoeff/SetDragCoeff) and by filters when updating
    /// estimated parameter values between iterations. Logs a warning (does
    /// nothing) if `key` is not a known parameter name.
    ///
    /// @param key    Parameter name to update (must already exist in `params_`).
    /// @param value  New scalar value for the parameter.
    void SetParam(const std::string& key, const Real value) {
      // set the parameter with the name=key to value
      std::vector<std::string> names = params_.GetNames();
      auto it = std::find(names.begin(), names.end(), key);
      if (it != names.end()) {
        int index = std::distance(names.begin(), it);
        // Defensive invariant check (optional but good)
        if (index >= static_cast<std::size_t>(params_.size())) {
          throw std::logic_error(
              fmt::format("[SetParam] Parameter index out of range for key '{}': "
                          "index = {}, params_.size() = {}, names.size() = {}. "
                          "This indicates a mismatch between parameter names and values.",
                          key, index, params_.size(), names.size()));
        }
        params_(index) = value;  // assuming value is a single element vector
      } else {
        spdlog::warn("Parameter {} not found in Dynamics parameters.", key);
      }
    }

    /// @brief Return the dynamics' current parameter vector (names + values).
    ParamState GetParams() const { return params_; }

    /// @brief Look up a single named parameter's value.
    ///
    /// @param key  Parameter name to look up.
    /// @return     The parameter's value, or `Real(0)` if `key` is not found.
    Real GetParam(const std::string& key) const {
      // get the parameter with the name=key
      std::vector<std::string> names = params_.GetNames();
      auto it = std::find(names.begin(), names.end(), key);
      if (it != names.end()) {
        int index = std::distance(names.begin(), it);
        return params_(index);
      }
      return Real(0);  // return zero if not found
    }

    /// @brief Return the `StateType` tag of the state vector this dynamics model
    /// operates on (e.g. `Cart6::TYPE`, `Attitude::TYPE`, `ClockState3::TYPE`).
    ///
    /// Used by state-converters and asset wiring code to verify that a dynamics
    /// object is paired with a compatible state representation.
    virtual StateType GetStateType() const = 0;
  };

}  // namespace lupnt
