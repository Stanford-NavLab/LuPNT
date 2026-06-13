/**
 * @file joint_orbit_clock_dynamics.h
 * @brief Coupled Cartesian orbit and clock dynamics.
 */

#pragma once

#include "lupnt/dynamics/clock_dynamics.h"
#include "lupnt/dynamics/numerical_orbit_dynamics.h"
#include "lupnt/environment/body.h"
#include "lupnt/states/state.h"

namespace lupnt {

  /**
   * @brief Numerical dynamics for a Cartesian orbit state coupled to a clock state.
   *
   * The state layout is [r, v, b, d] or [r, v, b, d, dr]. The orbit part is
   * propagated by an owned NumericalDynamics model. The clock part follows the
   * same bias-unit convention as ClockDynamics, with an optional relativistic
   * rate correction from the selected central body.
   */
  class JointOrbitClockDynamics : public NumericalDynamics {
  private:
    static constexpr int ORBIT_STATE_SIZE = 6;

    Ptr<NumericalDynamics> orbit_dynamics_;
    Ptr<ClockDynamics> clock_dynamics_;
    bool use_clock_relativity_ = true;
    bool relativity_center_body_set_ = false;
    BodyId relativity_center_body_ = BodyId::EARTH;
    Real reference_rate_offset_ = 0.0;
    bool add_clock_noise_ = false;
    Frame frame_ = Frame::UNDEFINED;
    UnitSystem units_ = SI_UNITS;

    /// @brief Pick the central body relative to which the relativistic clock-rate
    /// correction (RelativisticRate) is evaluated at time `t`.
    ///
    /// If `relativity_center_body_set_` (SetRelativityCenterBody), returns that
    /// fixed body. Else, if the integration frame's natural center is a real body
    /// (not SSB/Sun), returns that. Otherwise, for `NBodyDynamics` orbit models,
    /// returns whichever configured body the orbit state `orbit_state` is currently
    /// closest to; defaults to Earth if the orbit dynamics is not an
    /// `NBodyDynamics`.
    ///
    /// @param t            Time argument [s] (as passed to ComputeRates), used to
    ///                      evaluate candidate body positions via GetCenterState.
    /// @param orbit_state   Current orbital Cartesian state (position in
    ///                      `orbit_state.head(3)`), in the integration frame/units.
    /// @return             BodyId of the selected relativity center body.
    BodyId SelectRelativityCenter(Real t, const State& orbit_state) const;

    /// @brief Get the position/velocity of `center_body` relative to the
    /// integration frame's origin.
    ///
    /// Returns zero if `center_body` is already the frame's natural center;
    /// otherwise queries `GetBodyPosVel` at `t + GetLupntEpoch()` in the
    /// configured `frame_`/`units_`. Used by SelectRelativityCenter and
    /// RelativisticRate to re-center the orbit state on the chosen relativity
    /// center body.
    ///
    /// @param t            Time argument [s] (as passed to ComputeRates).
    /// @param center_body  Body to query the position/velocity of.
    /// @return             6-element `[r; v]` state of `center_body` relative to
    ///                      the integration frame origin, in `units_`.
    Vec6 GetCenterState(Real t, BodyId center_body) const;

    /// @brief Compute the special-relativistic clock-rate correction for the
    /// current orbit state.
    ///
    /// Selects the relativity center body (SelectRelativityCenter), re-centers the
    /// orbit state on it (GetCenterState), and evaluates
    /// `ClockDynamics::RelativisticRateCorrection` using that body's GM and the
    /// configured `reference_rate_offset_`. Called by ComputeRates when
    /// `use_clock_relativity_` is enabled to drive the clock-drift-rate term of the
    /// joint state derivative.
    ///
    /// @param t            Time argument [s] (as passed to ComputeRates).
    /// @param orbit_state   Current orbital Cartesian state in the integration
    ///                      frame/units.
    /// @return             Fractional relativistic clock-rate correction
    ///                      [dimensionless, i.e. ds/dt].
    Real RelativisticRate(Real t, const State& orbit_state) const;

    /// @brief Assemble a `JointOrbitClockState` from a propagated value vector,
    /// preserving/updating units and frame.
    ///
    /// Copies `values` into the head of a `JointOrbitClockState` derived from
    /// `x0`'s size/units, replaces the clock-portion units according to
    /// `clock_dynamics_->GetClockBiasUnit()`, and sets the frame to `frame_` if
    /// configured.
    ///
    /// @param x0      Reference state used to determine state size/initial units.
    /// @param values  New state values (orbit + clock portions).
    /// @return        `JointOrbitClockState` with `values`, appropriate units, and
    ///                frame.
    State BuildState(const State& x0, const VecX& values) const;

  public:
    /// @brief Construct with default (unset) orbit dynamics and a default
    /// `ClockDynamics`.
    ///
    /// Installs `ComputeRates` as the integration ODE via `SetODE`.
    JointOrbitClockDynamics();

    /// @brief Construct from a YAML configuration node.
    ///
    /// Reads the nested `orbit_dynamics` config (constructed via the
    /// `Dynamics`/`NumericalDynamics` asset factory and required to derive from
    /// `NumericalDynamics`), the nested `clock_dynamics` config, and the
    /// `clock_bias_unit`/`bias_unit`, `use_clock_relativity`/`use_relativity`,
    /// `relativity_center_body`, `reference_rate_offset`, `add_clock_noise`, and
    /// `frame` fields.
    explicit JointOrbitClockDynamics(Config& config);

    /// @brief Set the underlying orbit (position/velocity) dynamics model.
    ///
    /// Adopts the orbit dynamics' parameters (`SetParams`) and, if it is an
    /// `NBodyDynamics`, also adopts its integration frame and unit system.
    ///
    /// @param orbit_dynamics  Orbit dynamics model (must be non-null).
    void SetOrbitDynamics(Ptr<NumericalDynamics> orbit_dynamics);
    /// @brief Return the configured orbit dynamics model.
    Ptr<NumericalDynamics> GetOrbitDynamics() const { return orbit_dynamics_; }

    /// @brief Set the clock-bias propagation model (must be non-null).
    void SetClockDynamics(Ptr<ClockDynamics> clock_dynamics);
    /// @brief Return the configured clock-bias propagation model.
    Ptr<ClockDynamics> GetClockDynamics() const { return clock_dynamics_; }

    /// @brief Enable/disable the relativistic clock-rate coupling
    /// (RelativisticRate) in ComputeRates.
    void SetUseClockRelativity(bool use_clock_relativity) {
      use_clock_relativity_ = use_clock_relativity;
    }
    /// @brief Return whether relativistic clock-rate coupling is enabled.
    bool GetUseClockRelativity() const { return use_clock_relativity_; }

    /// @brief Fix the relativity center body used by SelectRelativityCenter,
    /// overriding automatic selection.
    void SetRelativityCenterBody(BodyId center_body) {
      relativity_center_body_ = center_body;
      relativity_center_body_set_ = true;
    }
    /// @brief Return the configured relativity center body (only meaningful if
    /// HasRelativityCenterBody() is true).
    BodyId GetRelativityCenterBody() const { return relativity_center_body_; }
    /// @brief Return whether a fixed relativity center body has been set via
    /// SetRelativityCenterBody.
    bool HasRelativityCenterBody() const { return relativity_center_body_set_; }
    /// @brief Clear a fixed relativity center body, reverting to automatic
    /// selection (SelectRelativityCenter).
    void ClearRelativityCenterBody() { relativity_center_body_set_ = false; }

    /// @brief Set a constant offset added to the relativistic clock-rate
    /// correction (see ClockDynamics::RelativisticRateCorrection).
    void SetReferenceRateOffset(Real reference_rate_offset) {
      reference_rate_offset_ = reference_rate_offset;
    }
    /// @brief Return the configured relativistic clock-rate reference offset.
    Real GetReferenceRateOffset() const { return reference_rate_offset_; }

    /// @brief Enable/disable adding stochastic clock process noise during
    /// Propagate.
    void SetAddClockNoise(bool add_clock_noise) { add_clock_noise_ = add_clock_noise; }
    /// @brief Return whether stochastic clock process noise is added during
    /// Propagate.
    bool GetAddClockNoise() const { return add_clock_noise_; }

    /// @brief Set the integration frame for the orbit portion of the joint state.
    void SetFrame(Frame frame) { frame_ = frame; }
    /// @brief Return the configured integration frame.
    Frame GetFrame() const { return frame_; }

    /// @brief Set the unit system used for the orbit portion of the joint state.
    void SetUnits(const UnitSystem& units) { units_ = units; }
    /// @brief Return the configured unit system.
    UnitSystem GetUnits() const { return units_; }

    using NumericalDynamics::Propagate;

    /// @brief Numerically propagate the coupled orbit+clock state over `[t0, tf]`.
    ///
    /// Integrates `ComputeRates` via `NumericalDynamics::Propagate` (skipping
    /// integration if `|tf - t0| < EPS`), then if `add_clock_noise_` and a clock
    /// model is configured, adds a sample of clock process noise
    /// (`clock_dynamics_->GetProcessNoise`) to the clock portion of the propagated
    /// state. Used by the filter time-update step (filters/) and by
    /// applications/ to propagate a satellite's coupled position/velocity/clock-bias
    /// state.
    ///
    /// @param x0  Initial `JointOrbitClockState` (`[r, v, b, d]` or
    ///            `[r, v, b, d, dr]`) at time `t0`.
    /// @param t0  Initial epoch [s, TDB since J2000].
    /// @param tf  Final epoch [s, TDB since J2000].
    /// @param u   Unused; must be nullptr (no control input supported).
    /// @return    Propagated `JointOrbitClockState` at time `tf`.
    State Propagate(const State& x0, Real t0, Real tf, const State* u = nullptr) override;

    using NumericalDynamics::PropagateWithParams;

    /// @brief Propagate the joint orbit+clock state after installing `params` on
    /// both this dynamics object and the underlying orbit dynamics.
    ///
    /// Calls `orbit_dynamics_->SetParams(params)` and `SetParams(params)`, then
    /// `Propagate(x0, t0, tf, u)`. Used by filters estimating orbit force-model
    /// parameters (e.g. SRP/drag coefficients) jointly with the orbit+clock state.
    ///
    /// @param x0      Initial `JointOrbitClockState` at time `t0`.
    /// @param t0      Initial epoch [s, TDB since J2000].
    /// @param tf      Final epoch [s, TDB since J2000].
    /// @param params  Orbit-dynamics parameter values/names to install before
    ///                 propagating.
    /// @param u       Unused; must be nullptr.
    /// @return        Propagated `JointOrbitClockState` at time `tf`.
    State PropagateWithParams(const State& x0, Real t0, Real tf, const ParamState& params,
                              const State* u = nullptr) override;

    /// @brief Compute the time derivative of the coupled orbit+clock state.
    ///
    /// Splits `x` into the orbit portion (`head(ORBIT_STATE_SIZE)`, tagged with
    /// `frame_` or `x`'s frame) and the clock portion, computes orbit rates via
    /// `orbit_dynamics_->ComputeRates`, computes the relativistic clock-rate term
    /// via RelativisticRate (if `use_clock_relativity_`), and forms the clock-state
    /// derivative as `[drift + relativistic_rate, drift_rate, 0]` (2-state:
    /// `[drift + relativistic_rate, 0]`). Installed as the integration ODE (via
    /// `SetODE` in the constructors) and called once per integration step by the
    /// integrators in numerics/integrator.h.
    ///
    /// @param t  Time argument [s] (as used by `ComputeRates`/`SetODE`, i.e.
    ///           relative to `GetLupntEpoch()`).
    /// @param x  Current `JointOrbitClockState` (8- or 9-element).
    /// @return   Time derivative of `x`, same layout/units as `x`.
    VecX ComputeRates(Real t, const State& x) const override;

    /// @brief Return `JointOrbitClockState::TYPE`, the coupled orbit+clock state
    /// type propagated by this model.
    StateType GetStateType() const override { return JointOrbitClockState::TYPE; }
  };

}  // namespace lupnt
