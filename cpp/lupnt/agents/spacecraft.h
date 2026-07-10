#pragma once

#include "lupnt/agents/agent.h"
#include "lupnt/dynamics/joint_orbit_clock_dynamics.h"

namespace lupnt {

  /// @brief A physical spacecraft whose *truth* is a joint orbit + clock 8-state
  /// `[r(3), v(3), clock_bias, clock_drift]`. It self-propagates that truth so
  /// other agents / sensors can query both its geometry and its clock (needed for
  /// two-way time/frequency transfer and pseudorange).
  ///
  /// Truth dynamics come entirely from the agent config, keeping the "onboard
  /// software never knows the truth" separation:
  ///  - **orbit**: the agent's `dynamics:` block (an `NBodyDynamics`-style force
  ///    model). When omitted, `Simulation::Setup` injects the shared
  ///    `world.force_model` here, so common orbit physics is declared once under
  ///    `world:`. Built with the state-transition matrix disabled (truth needs no STM).
  ///  - **clock**: an optional `clock:` block (`model`, `seed`, `add_noise`).
  ///
  /// The initial truth state is read from `initial_state:`, which accepts either an
  /// explicit `{r0_m, v0_mps}` (MOON_CI) or a `{class: ClassicalOE, ...}` block,
  /// plus optional `clock_bias_s` / `clock_drift_sps`.
  ///
  /// This is a general orbit+clock spacecraft agent, used by both the ISL
  /// scenarios and the Lunar-GNSS receiver. Truth is propagated in
  /// simulation-relative time; `JointOrbitClockDynamics` adds `GetLupntEpoch()`.
  class Spacecraft : public Agent {
  public:
    Spacecraft() = default;
    explicit Spacecraft(Config& config);

    void Setup() override { Agent::Setup(); }

    /// @brief Propagate the truth 8-state to `t` [s, sim-relative] (mutating) and log.
    void Step(Real t) override;

    /// @brief The truth joint 8-state `[r,v,clock_bias,clock_drift]` at time `t`
    /// (sim-relative [s]), without mutating the stored truth.
    State GetTruthStateAt(Real t) const;
    const State& GetTruthState() const { return truth_; }
    Real GetTruthTime() const { return truth_time_; }
    JointOrbitClockDynamics* GetTruthDynamics() const { return truth_dyn_.get(); }

    /// @brief Geometry facade: Cartesian `[r; v]` truth at `t`, in `Frame::MOON_CI`.
    Cart6 GetStateAt(Real t) const override;

  protected:
    /// @brief Build `truth_dyn_` (orbit from `dynamics:`, clock from `clock:`) and
    /// the initial `truth_` state from `initial_state:`. Shared by subclasses.
    void ConfigureTruth(Config& config);

    State truth_;
    Real truth_time_ = 0.0;
    Ptr<JointOrbitClockDynamics> truth_dyn_;
  };

}  // namespace lupnt
