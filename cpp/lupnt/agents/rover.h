#pragma once

#include "lupnt/agents/agent.h"

namespace lupnt {

  /// @brief A lunar-surface rover carrying an onboard navigation application.
  ///
  /// `Rover` is a thin `AgentWithDynamics` whose truth Cartesian state (and time) are set
  /// externally by the surface-navigation application (`SurfaceRoverNavApp`) at each epoch --
  /// the kinematic driving path is generated inside the app, not integrated from `dynamics_`.
  /// Its role in the agent framework is to **host that navigation application** (attached via
  /// `Agent::SetApplication` / the `application:` config block); the app reaches its owning
  /// rover through `Application::GetAgent()` and the shared environment through
  /// `Agent::GetWorld()`. `GetStateAt` returns the stored truth state.
  class Rover : public AgentWithDynamics {
  public:
    /// @brief Default-construct a bare `Rover` (no dynamics/state); used to host a
    /// `SurfaceRoverNavApp` estimator on the agent (e.g. the legacy `RunSurfaceNav` driver).
    Rover() = default;

    /// @brief Construct from a `rover` agent config block: reads `name`/`frequency`, then
    /// creates and attaches the `application:` (e.g. `SurfaceRoverNavApp`). No `dynamics:` block
    /// is required -- the hosted app drives the truth state.
    explicit Rover(Config& config);

    /// @brief Return the rover's stored truth Cartesian state `[r; v]` (set each epoch by the
    /// hosted app via `SetState`); the query time is ignored (no dynamics are integrated).
    Cart6 GetStateAt(Real /*t*/) const override { return state_; }
  };

}  // namespace lupnt
