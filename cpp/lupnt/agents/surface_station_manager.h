#pragma once

#include "lupnt/agents/agent.h"

namespace lupnt {

  /// @brief Ground-segment estimator agent for the distributed ISL ODTS scenario.
  ///
  /// A non-physical agent representing the surface-station network's central processing node.
  /// It hosts a `GroundOdtsApp` (its `application:` block) that runs the *centralized* ground
  /// filter -- a single EKF estimating every satellite's orbit + clock from the surface
  /// stations' one-way pseudoranges alone (no inter-satellite links), as a baseline against the
  /// satellites' distributed onboard filters. Having no state of its own, its `GetStateAt`
  /// returns a zero state (mirroring `GroundStationManager`/`IslOdtsManager`).
  class SurfaceStationManager : public Agent {
  public:
    SurfaceStationManager() = default;
    explicit SurfaceStationManager(Config& config) : Agent(config) {}

    Cart6 GetStateAt(Real /*t*/) const override { return Cart6(Vec6::Zero(), Frame::MOON_CI); }
  };

}  // namespace lupnt
