#pragma once

#include "lupnt/agents/agent.h"
#include "lupnt/dynamics/joint_orbit_clock_dynamics.h"

namespace lupnt {

  /// @brief A relay/navigation satellite whose *truth* is a joint orbit + clock 8-state
  /// `[r(3), v(3), clock_bias, clock_drift]`, for the distributed ISL ODTS scenario.
  ///
  /// Unlike the plain `Satellite` (a `Cart6` orbit-only agent), `IslSatellite` self-propagates
  /// a `JointOrbitClockDynamics` truth (orbit gravity field + stochastic clock), so neighbours
  /// can query both its geometry and its *clock* truth (needed for two-way time/frequency
  /// transfer). It hosts its onboard estimator as its `Application` (a `SatelliteOdtsApp`), which
  /// pulls its own and its neighbours' truth via `GetJointStateAt`. Truth is propagated in
  /// simulation-relative time; `JointOrbitClockDynamics` adds `GetLupntEpoch()` internally.
  class IslSatellite : public Agent {
  public:
    IslSatellite() = default;
    explicit IslSatellite(Config& config);

    /// @brief Schedule the periodic truth `Step` and set up the hosted application.
    void Setup() override { Agent::Setup(); }

    /// @brief Propagate the truth 8-state to `t` [s, sim-relative] (mutating) and log.
    void Step(Real t) override;

    /// @brief The truth joint 8-state `[r,v,clock_bias,clock_drift]` at time `t`
    /// (sim-relative [s]), without mutating the stored truth.
    State GetTruthStateAt(Real t) const;
    /// @brief The current stored truth 8-state and its time.
    const State& GetTruthState() const { return truth_; }
    Real GetTruthTime() const { return truth_time_; }
    JointOrbitClockDynamics* GetTruthDynamics() const { return truth_dyn_.get(); }

    /// @brief Geometry facade: Cartesian `[r; v]` truth at `t`, in `Frame::MOON_CI`.
    Cart6 GetStateAt(Real t) const override;

  private:
    State truth_;
    Real truth_time_ = 0.0;
    Ptr<JointOrbitClockDynamics> truth_dyn_;
  };

}  // namespace lupnt
