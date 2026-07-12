#pragma once

#include "lupnt/agents/agent.h"

namespace lupnt {

  /// @brief Coordinator agent for the ephemeris/almanac datasize-accuracy study (Example 9).
  ///
  /// A non-physical agent representing the navigation-message generation / ground processing
  /// node. It hosts an `EphemerisApp` (its `application:` block) that owns the whole study:
  /// propagating a lunar-satellite truth trajectory, fitting the `LansEphemeris` /
  /// `LansAlmanac` broadcast models over a sweep of fitting-window lengths, and sizing the
  /// broadcast bit budget. Having no state of its own, its `GetStateAt` returns a zero state
  /// (mirroring `SurfaceStationManager`).
  class EphemerisManager : public Agent {
  public:
    EphemerisManager() = default;
    explicit EphemerisManager(Config& config) : Agent(config) {}

    Cart6 GetStateAt(Real /*t*/) const override { return Cart6(Vec6::Zero(), Frame::MOON_CI); }
  };

}  // namespace lupnt
