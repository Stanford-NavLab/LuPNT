#pragma once
#include "lupnt/agents/agent.h"

namespace lupnt {

  class Simulation;

  class Constellation {
  protected:
    Simulation* sim_;
    std::vector<Ptr<Agent>> satellites_;
    std::string name_;

  public:
    /// @brief Construct an empty constellation (no satellites).
    Constellation();

    /// @brief Build a constellation of `Satellite` agents from a TLE file.
    ///
    /// Called by `AssetFactory<...>`-driven scenario setup to populate the
    /// `agents:` list from a `tle:` file: for each TLE entry, converts the
    /// mean elements to a `ClassicalOE` state (propagated to the current
    /// LuPNT epoch via the TLE's mean motion), builds a `CartesianTwoBodyDynamics`
    /// + `FixedPointingDynamics` `Satellite` with a GNSS transmitter device,
    /// and appends it to `satellites_`.
    /// @param agent_config YAML config node with `name`, `tle`, `frequency`,
    ///        and `attitude_dynamics` keys
    Constellation(Config& agent_config);

    /// @brief Call `Setup` on every satellite in the constellation. Invoked
    /// once by the simulation setup pass alongside individual agents'
    /// `Setup`.
    void Setup();

    /// @brief Call `Step` on every satellite in the constellation, advancing
    /// each to `time`. Invoked by the simulation's per-step scheduling.
    /// @param time Target simulation time [s]
    void Step(Real time);

    /// @brief Register this constellation (and all its satellites) with a
    /// `Simulation`, so each satellite can access the event scheduler and
    /// Cesium viewer via `Agent::GetSimulation`.
    void SetSimulation(Simulation* sim);
  };

}  // namespace lupnt
