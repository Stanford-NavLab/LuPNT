#pragma once

#include "lupnt/core/definitions.h"
#include "lupnt/devices/device.h"
#include "lupnt/dynamics/clock_dynamics.h"
#include "lupnt/states/state.h"

namespace lupnt {

  enum class ClockModel { OCXO, USO, CSAC, MINI_RAFS, RAFS, DSAC, UNDEFINED };

  /// @brief Onboard clock device: tracks a clock-bias state (bias, drift,
  /// drift-rate) using a `ClockDynamics` model and exposes the agent's
  /// current "true" time via `Read`.
  ///
  /// Attached to a satellite/rover/ground-station `Agent` (e.g. as the
  /// `"clock"` device) and used by `ClockDynamics`/`JointOrbitClockDynamics`
  /// to propagate clock error states for measurement simulation (pseudorange/
  /// carrier-phase clock-bias terms in `measurements/`) and by navigation
  /// filters as the truth/reference against which receiver clock states are
  /// estimated.
  class Clock : public Device {
  private:
    Real time_ = 0.0;  // [s]
    ClockState3 state_;
    ClockModel model_ = ClockModel::UNDEFINED;
    ClockBiasUnit bias_unit_ = ClockBiasUnit::SECONDS;
    bool use_relativity_ = false;
    bool relativity_center_body_set_ = false;
    BodyId relativity_center_body_ = BodyId::EARTH;
    UnitSystem relativity_units_ = SI_UNITS;
    Ptr<ClockDynamics> dynamics_;

    /// @brief Refresh the units (`ClockBiasUnit`-dependent strings) attached
    /// to `state_` after `bias_unit_` or the state size changes, so
    /// `GetState().GetUnits()` and `Log` report consistent unit labels.
    void ApplyClockStateUnits();

    /// @brief Build the `ClockRelativityContext` (position/velocity at `t0`
    /// and `tf` relative to the relativity center body, its GM, and the
    /// speed of light) needed by `ClockDynamics::PropagateWithRelativity`.
    ///
    /// Called from `Read` when `use_relativity_` is enabled. Queries the
    /// owning `Agent::GetStateAt` for the agent's Cartesian state at `t0` and
    /// `tf`, converts both to the inertial frame centered on
    /// `relativity_center_body_` (or the frame center of the state at `tf` if
    /// no center body was explicitly set), and packages the result together
    /// with that body's GM and `C` for the relativistic clock-rate
    /// correction.
    ///
    /// @param t0  Start time of the propagation interval [s]
    /// @param tf  End time of the propagation interval [s]
    /// @return    Context with centered Cartesian states, central-body GM
    ///            [m^3/s^2] (or consistent units), and speed of light `c`
    ClockRelativityContext BuildRelativityContext(Real t0, Real tf) const;

  public:
    Clock();
    Clock(Config& config);

    /// @brief Propagate the clock state from its last-read time to `t` and
    /// return the resulting clock time (true time `t` plus the current clock
    /// bias, converted from `bias_unit_` to seconds) [s].
    ///
    /// Called whenever a device/measurement model needs the agent's clock
    /// reading instead of true simulation time (e.g. when generating
    /// timestamped GNSS measurements). Internally advances `state_` via
    /// `ClockDynamics::Propagate` (or `PropagateWithRelativity` if
    /// `use_relativity_` is enabled, using `BuildRelativityContext`) from the
    /// previously stored `time_` to `t`, updates `time_ = t`, and adds the
    /// propagated bias (converted to seconds via
    /// `ClockDynamics::BiasUnitsToSeconds`) to `t`.
    ///
    /// @param t  Simulation (true) time to read the clock at [s]; must be
    ///           `>= ` the time of the previous `Read`/current `time_`
    /// @return   Clock-reading time = `t + bias` [s]
    Real Read(Real t);

    /// @brief Clock device step; currently a no-op (clock propagation happens
    /// lazily in `Read`).
    void Step(Real t) override;
    /// @brief Clock device setup; currently a no-op.
    void Setup() override;

    /// @brief Set the clock noise/drift model (e.g. `USO`, `CSAC`, `RAFS`)
    /// and propagate it to the underlying `ClockDynamics`.
    void SetModel(ClockModel model);
    /// @brief Get the configured clock noise/drift model.
    ClockModel GetModel() const { return model_; }

    /// @brief Set the clock-bias state (bias, drift, drift-rate) and refresh
    /// its unit labels via `ApplyClockStateUnits`.
    void SetState(const ClockState3& state) {
      state_ = state;
      ApplyClockStateUnits();
    }
    /// @brief Get the current clock-bias state (bias, drift, drift-rate), in
    /// units given by `GetClockBiasUnit`.
    ClockState3 GetState() const { return state_; }

    /// @brief Set the clock's internal last-propagated time [s] (does not
    /// propagate the state).
    void SetTime(Real time) { time_ = time; }
    /// @brief Get the clock's internal last-propagated time [s] (the time
    /// argument of the most recent `Read` call).
    Real GetTime() const { return time_; }

    /// @brief Set the unit (seconds, range-equivalent meters, or kilometers)
    /// in which the clock-bias state is represented, and propagate it to the
    /// underlying `ClockDynamics`.
    void SetClockBiasUnit(ClockBiasUnit unit);
    /// @brief Get the unit in which the clock-bias state is represented.
    ClockBiasUnit GetClockBiasUnit() const { return bias_unit_; }

    /// @brief Enable/disable general-relativistic clock-rate corrections
    /// (gravitational frequency shift + special-relativistic time dilation)
    /// in `Read`, via `ClockDynamics::PropagateWithRelativity`.
    void SetUseRelativity(bool use_relativity) { use_relativity_ = use_relativity; }
    /// @brief Whether relativistic clock-rate corrections are applied in `Read`.
    bool GetUseRelativity() const { return use_relativity_; }

    /// @brief Set the central body (e.g. `BodyId::EARTH`) whose
    /// gravitational potential and inertial frame define the relativistic
    /// clock-rate correction computed in `BuildRelativityContext`.
    void SetRelativityCenterBody(BodyId center_body) {
      relativity_center_body_ = center_body;
      relativity_center_body_set_ = true;
    }
    /// @brief Get the relativity center body (see `SetRelativityCenterBody`).
    BodyId GetRelativityCenterBody() const { return relativity_center_body_; }
    /// @brief Whether a relativity center body has been explicitly set; if
    /// false, `BuildRelativityContext` falls back to the frame center of the
    /// agent's state at the propagation end time.
    bool HasRelativityCenterBody() const { return relativity_center_body_set_; }
    /// @brief Clear a previously-set relativity center body (see
    /// `HasRelativityCenterBody`).
    void ClearRelativityCenterBody() { relativity_center_body_set_ = false; }

    /// @brief Set the unit system (e.g. SI) used to look up the relativity
    /// center body's physical constants (GM, speed of light) in
    /// `BuildRelativityContext`.
    void SetRelativityUnits(const UnitSystem& units) { relativity_units_ = units; }
    /// @brief Get the unit system used for relativity-related physical constants.
    UnitSystem GetRelativityUnits() const { return relativity_units_; }

    /// @brief Replace the underlying `ClockDynamics` model, propagating the
    /// current clock model and bias unit to it.
    void SetDynamics(Ptr<ClockDynamics> dynamics);
    /// @brief Get a raw (non-owning) pointer to the underlying `ClockDynamics`.
    ClockDynamics* GetDynamics() const { return dynamics_.get(); }
    /// @brief Get a shared pointer to the underlying `ClockDynamics`.
    Ptr<ClockDynamics> GetDynamicsSharedPtr() const { return dynamics_; }

    /// @brief Log the clock's current time and bias-state components
    /// (with unit-aware names, e.g. `"<name>/state/b_s"`) to the simulation's
    /// `DataLogger` output. Overrides `Device::Log`.
    virtual void Log(Real time) override;
  };

}  // namespace lupnt
