/**
 * @file clock.h
 * @author Stanford NAV LAB
 * @brief Clock class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <tuple>

#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"
#include "lupnt/dynamics/dynamics.h"
#include "lupnt/states/state.h"

namespace lupnt {

  enum class ClockModel;
  enum class ClockBiasUnit { SECONDS, METERS, KILOMETERS };

  struct ClockRelativityContext {
    Cart6 rv0_centered;
    Cart6 rvf_centered;
    Real GM = 0.0;
    Real c = C;
    Real reference_rate_offset = 0.0;
  };

  /// @brief Clock-bias/drift/drift-rate error-state propagation and noise model.
  ///
  /// `ClockDynamics` propagates a 2-state (`[bias, drift]`) or 3-state
  /// (`[bias, drift, drift-rate]`) `ClockState3`-type state using a polynomial state
  /// transition matrix (TwoStatePhi/ThreeStatePhi) plus, if a `ClockModel` oscillator
  /// model is set, additive random-walk process noise (TwoStateNoise/ThreeStateNoise)
  /// drawn via `GetProcessNoise`. Optionally applies a special-relativistic clock
  /// rate correction (PropagateWithRelativity / RelativisticRateCorrection). Used by
  /// `devices/clock.h::Clock` to propagate each agent's "true" clock-error state,
  /// and by `JointOrbitClockDynamics` to propagate the clock portion of the coupled
  /// orbit+clock filter state.
  class ClockDynamics : public Dynamics {
  protected:
    int seed_ = 0;
    Ptr<std::mt19937> rng_ = nullptr;
    ClockModel model_;
    bool add_noise_ = true;
    ClockBiasUnit bias_unit_ = ClockBiasUnit::SECONDS;

  private:
    /// @brief Shared implementation for the `Propagate`/`PropagateWithRelativity`
    /// overloads.
    ///
    /// Applies the bias-unit state transition matrix (TwoStatePhi/ThreeStatePhi),
    /// optionally adds oscillator process noise (if `model_ != UNDEFINED` and
    /// `add_noise_`), optionally applies the relativistic rate correction from
    /// `relativity` (averaging the rate at `t0` and `tf`), and updates the
    /// returned state's units via `GetStateUnits`.
    ///
    /// @param x0          Initial clock state (2- or 3-element `ClockState3`) at
    ///                     time `t0`, in `bias_unit_` units.
    /// @param t0          Initial epoch [s].
    /// @param tf          Final epoch [s].
    /// @param u           Unused; must be nullptr (clock dynamics has no control
    ///                     input).
    /// @param relativity  Optional relativistic-correction context (centered
    ///                     position/velocity at `t0` and `tf`, central-body GM,
    ///                     speed of light, reference rate offset); nullptr to skip
    ///                     the correction.
    /// @param stm         Output: state transition matrix (Phi_clk); left
    ///                     unmodified if nullptr.
    /// @return            Propagated clock state at time `tf`, in `bias_unit_`
    ///                     units.
    State PropagateImpl(const State& x0, Real t0, Real tf, const State* u,
                        const ClockRelativityContext* relativity, MatXd* stm);

  public:
    /// @brief Construct with an undefined (no-noise) clock model.
    ClockDynamics();

    /// @brief Construct from a YAML configuration node.
    ///
    /// Reads the optional `model` (ClockModel) and `clock_bias_unit`/`bias_unit`
    /// (ClockBiasUnit) fields.
    ClockDynamics(Config& config);

    /// @brief Look up the oscillator power-spectral-density coefficients (q1, q2,
    /// q3) for a given clock model.
    ///
    /// These are the white-frequency, random-walk-frequency, and (for 3-state
    /// models) frequency-drift PSD coefficients used by
    /// TwoStateNoise/ThreeStateNoise to build the process-noise covariance over a
    /// time step `dt`. Values are taken from manufacturer/mission specifications
    /// (e.g. OCXO, USO, CSAC, MINI_RAFS, RAFS, DSAC); see the source comments for
    /// references.
    ///
    /// @param clk_model  Clock oscillator model.
    /// @return           Tuple `(q1 [s^2/s], q2 [1/s], q3 [1/s^3])` (units as
    ///                    consumed by TwoStateNoise/ThreeStateNoise).
    static std::tuple<double, double, double> GetClockValues(ClockModel clk_model);

    /// @brief Return the multiplicative scale factor converting a clock bias from
    /// seconds to `unit`.
    ///
    /// `SECONDS -> 1`, `METERS -> C` (speed of light), `KILOMETERS -> C * KM_M`.
    /// Used to convert clock-bias states/covariances between time and
    /// range-equivalent units for pseudorange-based measurements
    /// (measurements/).
    ///
    /// @param unit  Target clock-bias unit.
    /// @return      Scale factor such that `value_in_unit = value_in_seconds * scale`.
    static double SecondsToBiasUnitScale(ClockBiasUnit unit);

    /// @brief Convert a clock-bias value from seconds to the given bias unit.
    ///
    /// @param value_s  Clock bias [s].
    /// @param unit     Target clock-bias unit.
    /// @return         Clock bias in `unit` (seconds, meters, or kilometers).
    static Real SecondsToBiasUnits(Real value_s, ClockBiasUnit unit);

    /// @brief Convert a clock-bias value from the given bias unit back to seconds.
    ///
    /// @param value  Clock bias in `unit`.
    /// @param unit   Source clock-bias unit.
    /// @return       Clock bias [s].
    static Real BiasUnitsToSeconds(Real value, ClockBiasUnit unit);

    /// @brief Return the unit strings for a clock state vector of the given size and
    /// bias unit.
    ///
    /// E.g. for `unit = SECONDS` and `state_size = 3`, returns
    /// `{"s", "s/s", "s/s^2"}`. Used to tag `ClockState3`/`State` objects with
    /// human-readable units after propagation.
    ///
    /// @param state_size  Clock state size (must be 2 or 3).
    /// @param unit        Clock-bias unit.
    /// @return            Vector of `state_size` unit strings (bias, drift,
    ///                    [drift-rate]).
    static std::vector<std::string> GetStateUnits(int state_size, ClockBiasUnit unit);

    /// @brief Compute the special/general-relativistic fractional clock-rate
    /// correction for an object at a given position/velocity relative to a central
    /// body.
    ///
    /// Returns `reference_rate_offset - (GM/r + 0.5*v^2) / c^2`, the leading-order
    /// gravitational + velocity time-dilation term (combined Sagnac/Shapiro-style
    /// rate offset) relative to a clock at infinity at rest. Used by
    /// PropagateImpl/PropagateWithRelativity to add a relativistic drift to the
    /// clock-bias rate, and by `JointOrbitClockDynamics::RelativisticRate` to couple
    /// the orbit state to the clock-bias rate in the joint filter state.
    ///
    /// @param r_centered           Position relative to the central body [m]
    ///                              (inertial frame, e.g. GCRF).
    /// @param v_centered           Velocity relative to the central body [m/s]
    ///                              (inertial frame, e.g. GCRF).
    /// @param GM                   Gravitational parameter of the central body
    ///                              [m^3/s^2].
    /// @param c                    Speed of light [m/s] (defaults to the global
    ///                              constant `C`).
    /// @param reference_rate_offset  Constant offset added to the result, e.g. to
    ///                              reference the rate to a different clock/frame
    ///                              [dimensionless].
    /// @return                     Fractional clock-rate correction [dimensionless,
    ///                              i.e. ds/dt].
    static Real RelativisticRateCorrection(const Vec3& r_centered, const Vec3& v_centered, Real GM,
                                           Real c = C, Real reference_rate_offset = 0.0);

    /// @brief Set the random-number-generator seed used for process-noise sampling.
    ///
    /// Reseeds the internal Mersenne Twister RNG (`rng_`).
    void SetSeed(int seed) {
      seed_ = seed;
      rng_ = MakePtr<std::mt19937>(seed_);
    }
    /// @brief Return the current RNG seed.
    int GetSeed() const { return seed_; }

    /// @brief Set the oscillator model used for process-noise generation.
    void SetModel(ClockModel model) { model_ = model; }
    /// @brief Return the configured oscillator model.
    ClockModel GetModel() const { return model_; }
    /// @brief Enable/disable adding stochastic process noise during propagation.
    void SetAddNoise(bool add_noise) { add_noise_ = add_noise; }
    /// @brief Return whether stochastic process noise is added during propagation.
    bool GetAddNoise() const { return add_noise_; }
    /// @brief Set the unit (seconds/meters/kilometers) used for the clock-bias
    /// state.
    void SetClockBiasUnit(ClockBiasUnit unit) { bias_unit_ = unit; }
    /// @brief Return the unit used for the clock-bias state.
    ClockBiasUnit GetClockBiasUnit() const { return bias_unit_; }

    /// @brief Draw a sample of the clock process noise over a time step `dt`.
    ///
    /// Builds the noise covariance via TwoStateNoise/ThreeStateNoise (using
    /// `model_` and `bias_unit_`) and draws one sample via `SampleMvNormal` using
    /// the internal RNG (`rng_`).
    ///
    /// @param dt          Time step [s].
    /// @param state_size  Clock state size (2 or 3).
    /// @return            Process-noise sample vector, in `bias_unit_` units.
    VecX GetProcessNoise(Real dt, int state_size);

    /// @brief Return the 2-state (bias, drift) clock state transition matrix for
    /// time step `dt`: `[[1, dt], [0, 1]]`.
    static Mat2 TwoStatePhi(Real dt);

    /// @brief Return the 3-state (bias, drift, drift-rate) clock state transition
    /// matrix for time step `dt`: `[[1, dt, dt^2/2], [0, 1, dt], [0, 0, 1]]`.
    static Mat3 ThreeStatePhi(Real dt);

    /// @brief Compute the 2-state clock process-noise covariance over time step
    /// `dt`, in seconds-equivalent units.
    ///
    /// Built from the white-frequency (`q1`) and random-walk-frequency (`q2`)
    /// power-spectral densities returned by GetClockValues for `clk_model`.
    ///
    /// @param clk_model  Clock oscillator model.
    /// @param dt         Time step [s].
    /// @return           2x2 process-noise covariance `Q` for `[bias, drift]`,
    ///                    in seconds^2 / (seconds/s)^2.
    static Mat2 TwoStateNoise(ClockModel clk_model, Real dt);

    /// @brief 2-state process-noise covariance, scaled to the given clock-bias
    /// unit.
    ///
    /// Equivalent to `TwoStateNoise(clk_model, dt) * scale^2` where `scale =
    /// SecondsToBiasUnitScale(unit)`.
    static Mat2 TwoStateNoise(ClockModel clk_model, Real dt, ClockBiasUnit unit);

    /// @brief Compute the 3-state clock process-noise covariance over time step
    /// `dt`, in seconds-equivalent units.
    ///
    /// Built from the white-frequency (`q1`), random-walk-frequency (`q2`), and
    /// frequency-drift (`q3`) power-spectral densities returned by GetClockValues
    /// for `clk_model`.
    ///
    /// @param clk_model  Clock oscillator model.
    /// @param dt         Time step [s].
    /// @return           3x3 process-noise covariance `Q` for
    ///                    `[bias, drift, drift-rate]`, in seconds-equivalent units.
    static Mat3 ThreeStateNoise(ClockModel clk_model, Real dt);

    /// @brief 3-state process-noise covariance, scaled to the given clock-bias
    /// unit.
    ///
    /// Equivalent to `ThreeStateNoise(clk_model, dt) * scale^2` where `scale =
    /// SecondsToBiasUnitScale(unit)`.
    static Mat3 ThreeStateNoise(ClockModel clk_model, Real dt, ClockBiasUnit unit);

    using Dynamics::Propagate;

    /// @brief Propagate the clock-bias state and compute its state transition
    /// matrix, without relativistic correction.
    ///
    /// Equivalent to `PropagateImpl(x0, t0, tf, u, nullptr, stm)`.
    ///
    /// @param x0   Initial clock state (2- or 3-element) at time `t0`.
    /// @param t0   Initial epoch [s].
    /// @param tf   Final epoch [s].
    /// @param u    Unused; must be nullptr.
    /// @param stm  Output: clock state transition matrix; left unmodified if
    ///             nullptr.
    /// @return     Propagated clock state at time `tf`.
    State Propagate(const State& x0, Real t0, Real tf, const State* u, MatXd* stm) override;

    /// @brief Propagate the clock-bias state without relativistic correction or
    /// STM output.
    ///
    /// Equivalent to `PropagateImpl(x0, t0, tf, u, nullptr, nullptr)`.
    State Propagate(const State& x0, Real t0, Real tf, const State* u = nullptr) override;

    /// @brief Propagate the clock-bias state including a special-relativistic clock
    /// rate correction.
    ///
    /// Called by `JointOrbitClockDynamics`/`devices/clock.h::Clock` when the
    /// object's orbital position/velocity relative to a central body is known, so
    /// that the gravitational + velocity time-dilation rate (see
    /// RelativisticRateCorrection) is added to the clock-bias drift over `[t0,
    /// tf]`.
    ///
    /// @param x0          Initial clock state at time `t0`.
    /// @param t0          Initial epoch [s].
    /// @param tf          Final epoch [s].
    /// @param relativity  Relativistic-correction context: centered
    ///                     position/velocity at `t0` and `tf`, central-body GM,
    ///                     speed of light, and reference rate offset.
    /// @param u           Unused; must be nullptr.
    /// @param stm         Output: clock state transition matrix; left unmodified if
    ///                     nullptr.
    /// @return            Propagated clock state at time `tf`, including the
    ///                     relativistic bias contribution.
    State PropagateWithRelativity(const State& x0, Real t0, Real tf,
                                  const ClockRelativityContext& relativity,
                                  const State* u = nullptr, MatXd* stm = nullptr);

    /// @brief Return `ClockState3::TYPE`, the clock-bias state type propagated by
    /// this model.
    StateType GetStateType() const override { return ClockState3::TYPE; }
  };
}  // namespace lupnt
