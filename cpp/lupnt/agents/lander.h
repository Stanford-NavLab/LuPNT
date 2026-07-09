#pragma once

#include "lupnt/agents/agent.h"

namespace lupnt {

  /// @brief A lunar lander performing powered descent to the surface.
  ///
  /// `Lander` is a thin `AgentWithDynamics` whose truth Cartesian state (and time) are set
  /// externally by the descent simulation (`RunLanderNav`) at each epoch. Its purpose in the
  /// agent framework is to **host the onboard navigation application**: a `LanderNavApp`
  /// (attached via `Agent::SetApplication`) runs the lander's error-state INS EKF, fusing the
  /// IMU, radar altimeter, crater bearings, and LunaNet pseudoranges. The app reaches its
  /// owning lander through `Application::GetAgent()`.
  ///
  /// The descent path is generated kinematically (not by integrating `dynamics_`), so
  /// `GetStateAt` simply returns the stored truth state rather than propagating a dynamics
  /// model.
  class Lander : public AgentWithDynamics {
  public:
    Lander() = default;
    explicit Lander(Config& agent_config) : AgentWithDynamics(agent_config) {}

    /// @brief Return the lander's stored truth Cartesian state `[r; v]`.
    /// The descent simulation sets this each epoch via `SetState`; no dynamics model is
    /// integrated, so the query time is ignored.
    Cart6 GetStateAt(Real /*t*/) const override { return state_; }
  };

}  // namespace lupnt
