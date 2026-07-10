#pragma once

#include "lupnt/agents/agent.h"

namespace lupnt {

  /// @brief Ground-segment estimator agent.
  ///
  /// A non-physical agent representing the mission's ground control / processing
  /// node. It hosts a `GroundStationManagerApp` (its `application:` block) that
  /// aggregates the observations reported by the `GroundStationTrackingApp`s on the
  /// individual `GroundStation` agents and runs the centralized orbit-determination
  /// filter. Having no state of its own, its `GetStateAt` returns a zero state.
  class GroundStationManager : public Agent {
  public:
    GroundStationManager() = default;
    explicit GroundStationManager(Config& config) : Agent(config) {}

    Cart6 GetStateAt(Real /*t*/) const override { return Cart6(Vec6::Zero(), Frame::MOON_CI); }
  };

}  // namespace lupnt
