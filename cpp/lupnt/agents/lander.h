#pragma once

#include "lupnt/agents/agent.h"

namespace lupnt {

  /// @brief A lunar lander performing powered descent to the surface.
  ///
  /// `Lander` is a thin `AgentWithDynamics` whose truth Cartesian state (and time) are set
  /// externally by the descent navigation application (`LanderNavApp`) at each epoch -- the
  /// descent path is generated inside the app, not integrated from `dynamics_`. Its role in the
  /// agent framework is to **host that navigation application** (attached via
  /// `Agent::SetApplication` / the `application:` config block); the app reaches its owning
  /// lander through `Application::GetAgent()` and the shared environment through
  /// `Agent::GetWorld()`. `GetStateAt` returns the stored truth state.
  class Lander : public AgentWithDynamics {
  public:
    Lander() = default;

    /// @brief Construct from a `lander` agent config block: reads `name`/`frequency`, then
    /// creates and attaches the `application:` (e.g. `LanderNavApp`). No `dynamics:` block is
    /// required -- the hosted app drives the truth state.
    explicit Lander(Config& config);

    /// @brief Return the lander's stored truth Cartesian state `[r; v]` (set each epoch by the
    /// hosted app via `SetState`); the query time is ignored (no dynamics are integrated).
    Cart6 GetStateAt(Real /*t*/) const override { return state_; }
  };

}  // namespace lupnt
