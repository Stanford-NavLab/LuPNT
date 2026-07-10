#pragma once

#include "lupnt/agents/agent.h"

namespace lupnt {

  /// @brief Coordinator agent for the lunar-orbiting GNSS ODTS scenario (Example 6).
  ///
  /// A non-physical agent representing the receiver's ground/processing node. It hosts a
  /// `LunarGnssOdtsApp` (its `application:` block) that owns the whole GNSS ODTS pipeline:
  /// the receiver truth trajectory, the cislunar GNSS sidelobe link geometry / CN0, and the
  /// UDU EKF (or UDU stochastic-cloning EKF when TDCP is enabled) Monte Carlo run. Having no
  /// state of its own, its `GetStateAt` returns a zero state (mirroring `IslOdtsManager`).
  class LunarGnssManager : public Agent {
  public:
    LunarGnssManager() = default;
    explicit LunarGnssManager(Config& config) : Agent(config) {}

    Cart6 GetStateAt(Real /*t*/) const override { return Cart6(Vec6::Zero(), Frame::MOON_CI); }
  };

}  // namespace lupnt
