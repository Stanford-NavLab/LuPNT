#pragma once

#include <string>
#include <vector>

#include "lupnt/agents/agent.h"
#include "lupnt/agents/spacecraft.h"

namespace lupnt {

  /// @brief A lunar navigation constellation as a SINGLE agent that stands in for N similar
  /// navigation satellites, so a scenario need not declare an individual agent per satellite.
  ///
  /// It builds N child `Spacecraft` (each self-propagating its orbit+clock truth) from one config
  /// entry, in either of two ways:
  ///  - **explicit** `satellites:` — a list of per-satellite `{name, initial_state}` (or the
  ///    `{name, r0_m, v0_mps}` shorthand); every satellite inherits the shared `dynamics:` (and
  ///    optional `clock:`) template, so the force model is written once.
  ///  - **Walker** `walker:` — a symmetric Walker constellation (`n_planes` x `sats_per_plane`)
  ///    of a common frozen orbit (`a`, `e`, `i`, `omega`) in `frame`, with the RAAN spread evenly
  ///    over the planes and the mean anomaly evenly within each plane (plus an inter-plane phase).
  ///
  /// Consumers (e.g. a surface rover / lander nav app) resolve this agent by name and query
  /// `GetSatelliteStateAt(j, t)` — they no longer propagate the relays themselves.
  class LunarNavConstellation : public Agent {
  public:
    LunarNavConstellation() = default;
    explicit LunarNavConstellation(Config& config);

    /// @brief Register each child with the owning simulation and set it up.
    void Setup() override;
    /// @brief Propagate/log each child satellite's truth.
    void Step(Real t) override;

    /// @brief Cartesian state `[r; v]` of the whole constellation is undefined; returns the first
    /// satellite's state (or a zero state if empty). Prefer `GetSatelliteStateAt`.
    Cart6 GetStateAt(Real t) const override;

    /// @brief Number of satellites in the constellation.
    int NumSatellites() const { return static_cast<int>(sats_.size()); }
    /// @brief Name of satellite `j`.
    const std::string& SatelliteName(int j) const { return sat_names_[j]; }
    /// @brief Cartesian truth state `[r; v]` (Frame::MOON_CI) of satellite `j` at time `t` [s].
    Cart6 GetSatelliteStateAt(int j, Real t) const { return sats_[j]->GetStateAt(t); }

  private:
    /// @brief Build one child Spacecraft from the shared `dynamics:`/`clock:` template plus the
    /// given `initial_state` YAML node and name.
    void AddSatellite(const Config& shared, const YAML::Node& initial_state, const std::string& nm);

    std::vector<Ptr<Spacecraft>> sats_;
    std::vector<std::string> sat_names_;
  };

}  // namespace lupnt
