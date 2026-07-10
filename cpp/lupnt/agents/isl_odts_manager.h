#pragma once

#include "lupnt/agents/agent.h"

namespace lupnt {

  /// @brief Coordinator agent for the distributed inter-satellite-link (ISL) ODTS scenario.
  ///
  /// A non-physical agent representing the constellation-level coordination / processing
  /// node. It hosts an `IslOdtsCoordinatorApp` (its `application:` block) that owns the whole
  /// per-epoch ODTS body: truth propagation, crosslink + surface-station measurement
  /// synthesis, the N parallel onboard Schmidt-EKF filters, the centralized ground filter,
  /// and the periodic consider-state exchange. Having no state of its own, its `GetStateAt`
  /// returns a zero state (mirroring `GroundStationManager`).
  class IslOdtsManager : public Agent {
  public:
    IslOdtsManager() = default;
    explicit IslOdtsManager(Config& config) : Agent(config) {}

    Cart6 GetStateAt(Real /*t*/) const override { return Cart6(Vec6::Zero(), Frame::MOON_CI); }
  };

}  // namespace lupnt
