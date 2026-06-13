/**
 * @file integrator.h
 * @author Stanford NAV LAB
 * @brief Integrator interfaces
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <functional>

#include "lupnt/core/constants.h"
#include "lupnt/core/error.h"
#include "lupnt/core/object.h"
#include "lupnt/states/state.h"

namespace lupnt {

  /// @brief Right-hand side of an ODE: `dx/dt = f(t, x)` for use with the `Integrator`
  /// hierarchy below.
  ///
  /// Numerical dynamics models (e.g. `NumericalOrbitDynamics`, `ClockDynamics`,
  /// `JointOrbitClockDynamics`) build an `ODE` from their force/derivative models and pass
  /// it to `Integrator::Propagate`/`PropagateEx` to advance the state.
  using ODE = std::function<VecX(Real, const State&)>;

  /// @brief Selects which concrete `Integrator` subclass `NumericalOrbitDynamics` (and
  /// similar dynamics classes) should instantiate via `SetIntegrator`.
  enum class IntegratorType {
    RK4,
    RK8,
    RKF45,
    PD45,
  };
  /// @brief Default integrator used when a dynamics model does not explicitly call
  /// `SetIntegrator`.
  constexpr IntegratorType default_integrator = IntegratorType::RK4;

  /// @brief Tolerances, iteration limits, and optional early-termination predicate shared
  /// by all `Integrator` subclasses.
  ///
  /// Configured by dynamics models via `Integrator::SetParams` (e.g.
  /// `NumericalOrbitDynamics::SetIntegratorParams`) to control adaptive step-size
  /// integrators (`IRKF`/`RKF45`, `PD45`) and to allow a propagation to stop early when a
  /// user-defined condition on the state is met.
  class IntegratorParams {
  public:
    int max_iter = 20;
    double abstol = 1e-6;
    double reltol = 1e-6;

    // User-specified termination: return true to stop
    std::function<bool(Real, const VecX&)> terminate_if = nullptr;

    IntegratorParams() = default;
    IntegratorParams(int max_iter, double abstol, double reltol)
        : max_iter(max_iter), abstol(abstol), reltol(reltol) {
      CheckIntegratorParams();
    };
    IntegratorParams(int max_iterm, double abstol, double reltol,
                     std::function<bool(Real, const VecX&)> terminate_if)
        : max_iter(max_iterm),
          abstol(abstol),
          reltol(reltol),
          terminate_if(std::move(terminate_if)) {
      CheckIntegratorParams();
    };

    /// @brief Validate that `max_iter`, `abstol`, and `reltol` are positive.
    ///
    /// Called by every `IntegratorParams` constructor; throws (via `LUPNT_CHECK`) if the
    /// parameters are invalid.
    void CheckIntegratorParams();
  };

  /// @brief Why an `Integrator::PropagateEx` call stopped: it reached the requested final
  /// time, or the user-supplied `IntegratorParams::terminate_if` predicate fired.
  enum class TerminationReason { ReachedTf, UserCondition };

  /// @brief Result of `Integrator::PropagateEx`: the final state/time reached, why
  /// propagation stopped, and how many integration steps were taken.
  struct IntegratorResult {
    VecX x;
    Real t;
    TerminationReason reason;
    int steps;
  };

  /// @brief Abstract base class for fixed- and adaptive-step ODE integrators used to
  /// propagate `State`/`VecX` dynamics in time.
  ///
  /// Concrete subclasses (`RK4`, `RK8`, `RKF45`, `PD45`) implement `Step`; numerical
  /// dynamics models (e.g. `NumericalOrbitDynamics`, `ClockDynamics`,
  /// `JointOrbitClockDynamics`) own an `Integrator` instance (selected via
  /// `SetIntegrator`/`IntegratorType`) and call `Propagate`/`PropagateEx` once per
  /// `Dynamics::Propagate` call to advance the state from `t0` to `tf`.
  class Integrator {
  private:
  protected:
    bool print_progress_ = false;
    IntegratorParams params_;

  public:
    virtual ~Integrator() {};

    /// @brief Propagate a `State` from `t0` to `tf` under ODE `odefunc`, taking
    /// fixed/limited steps of (at most) size `dt`.
    ///
    /// Called by numerical dynamics models (e.g.
    /// `NumericalOrbitDynamics::Propagate`/`ClockDynamics::Propagate`) to advance the
    /// propagated state by one `Dynamics::Propagate` interval.
    ///
    /// @param odefunc  Right-hand side `dx/dt = f(t, x)`
    /// @param t0       Initial time [s]
    /// @param tf       Final time [s] (must be >= `t0`; steps are clamped so the last
    ///                  step lands exactly on `tf`)
    /// @param x0       Initial state
    /// @param dt       Nominal step size [s] (must be positive)
    /// @return         Propagated state at `tf`
    State Propagate(const ODE& odefunc, Real t0, Real tf, const State& x0, Real dt);

    /// @brief Same as `Propagate(odefunc, t0, tf, x0, dt)`, additionally computing the
    /// state-transition (sensitivity) matrix `d(xf)/d(x0)` via parallel finite-difference
    /// Jacobians (`JacobianParallel`) if `J` is non-null.
    ///
    /// @param J  Output state-transition matrix (size n x n, where n = `x0.size()`); pass
    ///           `nullptr` to skip the Jacobian computation
    State Propagate(const ODE& odefunc, Real t0, Real tf, const State& x0, Real dt, MatXd* J);

    /// @brief `VecX` overload of `Propagate(odefunc, t0, tf, x0, dt)`; internally delegates
    /// to `PropagateEx` and returns just the final state vector.
    VecX Propagate(const ODE& odefunc, Real t0, Real tf, const VecX& x0, Real dt);

    /// @brief `VecX` overload of the Jacobian-computing `Propagate` above.
    VecX Propagate(const ODE& odefunc, Real t0, Real tf, const VecX& x0, Real dt, MatXd* J);

    /// @brief Propagate `x0` from `t0` to `tf` under `odefunc`, stepping by (the magnitude
    /// of) `dt` in the direction implied by `sign(tf - t0)`, and report how/when
    /// propagation stopped.
    ///
    /// Steps are clamped so the last step lands exactly on `tf`. If
    /// `params_.terminate_if(t, x)` is set and returns true (checked once before stepping
    /// and again after each step), propagation stops early with
    /// `TerminationReason::UserCondition`. Used by
    /// `NumericalOrbitDynamics::PropagateEx`/`PropagateWithSTM`-style callers that need the
    /// actual stop time/step count, not just the final state.
    ///
    /// @param odefunc  Right-hand side `dx/dt = f(t, x)`
    /// @param t0       Initial time [s]
    /// @param tf       Target final time [s] (may be less than `t0` for backward
    ///                  propagation)
    /// @param x0       Initial state vector
    /// @param dt       Step size magnitude [s] (must be positive; direction is inferred
    ///                  from `tf - t0`)
    /// @return         Final state, time, termination reason, and step count
    IntegratorResult PropagateEx(const ODE& odefunc, Real t0, Real tf, const VecX& x0, Real dt);

    /// @brief Same as `PropagateEx(odefunc, t0, tf, x0, dt)`, additionally computing the
    /// state-transition matrix `d(xf)/d(x0)` via `JacobianParallel` if `J` is non-null.
    ///
    /// @param J  Output state-transition matrix (size n x n, where n = `x0.size()`); pass
    ///           `nullptr` to skip the Jacobian computation
    IntegratorResult PropagateEx(const ODE& odefunc, Real t0, Real tf, const VecX& x0, Real dt,
                                 MatXd* J);

    /// @brief Enable/disable a progress bar (via `Logger::GetProgressBar`) during
    /// `Propagate`/`PropagateEx` calls -- useful for long-running propagations in
    /// interactive/scripted runs.
    void SetPrintProgress(bool print) { print_progress_ = print; }

    /// @brief Advance `x` by one integration step of size `dt` under ODE `f` at time `t`.
    ///
    /// Implemented by each concrete integrator (`RK4`, `RK8`, `IRKF`/`RKF45`, `PD45`) and
    /// called repeatedly by `Propagate`/`PropagateEx` to march the state from `t0` to `tf`.
    ///
    /// @param f   Right-hand side `dx/dt = f(t, x)`
    /// @param t   Current time [s]
    /// @param x   Current state
    /// @param dt  Step size [s] (signed; adaptive integrators may modify it in-place to
    ///            reflect the step actually taken)
    /// @return    State after stepping by `dt`
    virtual State Step(const ODE& f, Real t, const State& x, Real dt) = 0;

    /// @brief Set the tolerances/iteration-limit/termination predicate used by
    /// `Propagate`/`PropagateEx` and adaptive `Step` implementations.
    void SetParams(IntegratorParams params) { params_ = params; };

    /// @brief Install (or clear, with `nullptr`) an early-termination predicate evaluated
    /// during `PropagateEx`; propagation stops with `TerminationReason::UserCondition`
    /// once `pred(t, x)` returns true.
    void SetTerminateIf(std::function<bool(Real, const VecX&)> pred) {
      params_.terminate_if = std::move(pred);
    }
  };

  /// @brief Classical 4th-order Runge-Kutta (RK4) fixed-step integrator.
  ///
  /// `default_integrator` for `NumericalOrbitDynamics` and other numerical dynamics
  /// models; a simple, computationally cheap fixed-step method appropriate when `dt` is
  /// small relative to the dynamics' time scales.
  class RK4 : public Integrator {
  public:
    /// @brief One step of classical 4th-order Runge-Kutta integration: evaluates `f` at 4
    /// stages spanning `[t, t+dt]` and averages the resulting derivatives.
    ///
    /// @param f   Right-hand side `dx/dt = f(t, x)`
    /// @param t   Current time [s]
    /// @param x   Current state
    /// @param dt  Step size [s]
    /// @return    State after stepping by `dt`
    State Step(const ODE& f, Real t, const State& x, Real dt);
  };

  /// @brief 8th-order Runge-Kutta fixed-step integrator (10-stage, Cooper-Verner-type
  /// coefficients).
  ///
  /// Higher-accuracy alternative to `RK4` for `NumericalOrbitDynamics` and similar models
  /// when larger step sizes or tighter accuracy are needed at the cost of more force-model
  /// evaluations per step.
  class RK8 : public Integrator {
  public:
    /// @brief One step of 8th-order Runge-Kutta integration: evaluates `f` at 10 stages
    /// spanning `[t, t+dt]` and combines the resulting derivatives with fixed weights.
    ///
    /// @param f   Right-hand side `dx/dt = f(t, x)`
    /// @param t   Current time [s]
    /// @param x   Current state
    /// @param dt  Step size [s]
    /// @return    State after stepping by `dt`
    State Step(const ODE& f, Real t, const State& x, Real dt);
  };

  /// @brief Abstract base for embedded Runge-Kutta-Fehlberg-type integrators with adaptive
  /// step-size control.
  ///
  /// `Step` repeatedly calls the subclass-provided `Update` to compute a low- and
  /// high-order solution pair, checks the relative error via `ComputeRelError`, and shrinks
  /// `dt` and retries (up to `IntegratorParams::max_iter` times) until the error is within
  /// tolerance. `RKF45` is the concrete 4(5)-order instantiation used by
  /// `NumericalOrbitDynamics` when `IntegratorType::RKF45` is selected.
  class IRKF : public Integrator {
  private:
    int order_;

  public:
    /// @brief Construct an embedded RKF-type integrator of the given low-order accuracy
    /// `order` (used by `ComputeRelError`'s step-size-control exponent).
    IRKF(int order) : order_(order) {};

    /// @brief One adaptive step: repeatedly calls `Update` and `ComputeRelError`,
    /// shrinking/growing `dt` until the embedded low/high-order solutions agree to within
    /// `IntegratorParams::abstol`/`reltol`, or `max_iter` is exceeded (throws via
    /// `LUPNT_CHECK` in that case).
    ///
    /// @param f   Right-hand side `dx/dt = f(t, x)`
    /// @param t   Current time [s]
    /// @param x   Current state
    /// @param dt  Step size [s]; updated in-place by `ComputeRelError`
    /// @return    Low-order solution after stepping by (the final) `dt`
    State Step(const ODE& f, Real t, const State& x, Real dt) override;

    /// @brief Compute the relative error norm between the embedded low- and high-order
    /// solutions and adjust `dt` accordingly (PI-type step-size controller).
    ///
    /// Called once per iteration inside `Step` after `Update` produces `x_new_low` and
    /// `x_new_high`.
    ///
    /// @param x_new_low   Lower-order embedded solution
    /// @param x_new_high  Higher-order embedded solution
    /// @param dt          Step size [s]; rescaled in-place by the step-size controller
    ///                     (clamped to [0.5, 2.0] times its input value) for the next
    ///                     attempt
    /// @return            True if the error is within `IntegratorParams::abstol`/`reltol`
    ///                     for every component (step accepted); false otherwise (caller
    ///                     should retry with the rescaled `dt`)
    bool ComputeRelError(const State& x_new_low, const State& x_new_high, Real& dt);

    /// @brief Compute the embedded low- and high-order solution pair for one step of size
    /// `dt`, using the subclass's Butcher tableau.
    ///
    /// @param f          Right-hand side `dx/dt = f(t, x)`
    /// @param t          Current time [s]
    /// @param x          Current state
    /// @param dt         Step size [s]
    /// @param x_new_low  Output: lower-order solution at `t + dt`
    /// @param x_new_high Output: higher-order solution at `t + dt`
    virtual void Update(const ODE& f, Real t, const State& x, Real dt, State& x_new_low,
                        State& x_new_high)
        = 0;
    virtual ~IRKF() = default;
  };

  /// @brief Runge-Kutta-Fehlberg 4(5) adaptive-step integrator.
  ///
  /// Concrete `IRKF` instantiation selected by `NumericalOrbitDynamics::SetIntegrator` via
  /// `IntegratorType::RKF45`; suited to propagations where step size should automatically
  /// shrink/grow to meet `IntegratorParams::abstol`/`reltol`.
  class RKF45 : public IRKF {
  public:
    RKF45() : IRKF(4) {};

    /// @brief Compute the RKF45 embedded 4th/5th-order solution pair (6-stage Butcher
    /// tableau) for one step of size `dt`.
    ///
    /// @param f          Right-hand side `dx/dt = f(t, x)`
    /// @param t          Current time [s]
    /// @param x          Current state
    /// @param dt         Step size [s]
    /// @param x_new_low  Output: 4th-order solution at `t + dt`
    /// @param x_new_high Output: 5th-order solution at `t + dt`
    void Update(const ODE& f, Real t, const State& x, Real dt, State& x_new_low,
                State& x_new_high) override;
  };

  /// @brief Dormand-Prince 4(5) ("PD45") adaptive-step integrator with self-contained
  /// step-size control (does not rely on `IRKF`/`ComputeRelError`).
  ///
  /// Alternative adaptive integrator selectable via `IntegratorType::PD45` in
  /// `NumericalOrbitDynamics::SetIntegrator`.
  class PD45 : public Integrator {
  private:
    static const std::array<std::array<double, 6>, 7> A_;
    static const std::array<double, 7> b_;
    static const std::array<double, 7> b_star_;

  public:
    /// @brief One adaptive Dormand-Prince 4(5) step: computes the 7-stage embedded
    /// 4th/5th-order solution pair, and shrinks `dt` and retries (up to
    /// `IntegratorParams::max_iter` times) until the scaled error norm is <= 1, or throws
    /// if `dt` underflows `1e-3` or `max_iter` is exceeded.
    ///
    /// @param f   Right-hand side `dx/dt = f(t, x)`
    /// @param t   Current time [s]; advanced in-place by the accepted `dt` on success
    /// @param x   Current state
    /// @param dt  Step size [s]; shrunk in-place on rejected attempts
    /// @return    5th-order ("high") solution after stepping by the accepted `dt`
    State Step(const ODE& f, Real t, const State& x, Real dt) override;
  };

}  // namespace lupnt
