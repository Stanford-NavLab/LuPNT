#pragma once

#include <string>

#include "lupnt/core/config.h"
#include "lupnt/core/definitions.h"
#include "lupnt/core/object.h"
#include "lupnt/interfaces/lola_dem.h"
#include "lupnt/states/state.h"

namespace lupnt {

  class Simulation;
  class NBodyDynamics;

  /// @brief Shared, read-only physical environment for a `Simulation`.
  ///
  /// The `World` owns the pieces of a scenario that are common to every agent and
  /// expensive to re-specify per agent: the simulation epoch, the reference frame,
  /// and the environment models. Two flavors of environment are supported (either
  /// or both may be present in the `world:` block):
  ///  - **Orbital**: a `force_model:` block (central bodies, gravity degree/order,
  ///    third bodies, SRP/relativity, integrator) → `MakeDynamics()` builds a
  ///    matching `NBodyDynamics`, so truth and estimator share one force model.
  ///  - **Surface**: a central-body point-mass `Gravity()` plus, optionally, a
  ///    `dem:` terrain block → a loaded `LunarDem` and its local East-North-Up
  ///    tangent frame, for surface-navigation scenarios (rover/lander).
  ///
  /// The `World` deliberately does **not** propagate agents or own their states:
  /// each `AgentWithDynamics` self-propagates via its own dynamics pointer. The
  /// `World` only *provides* the environment and a read-only truth facade
  /// (`GetStateAt`) so applications can ask "where is agent X at time t" without
  /// reaching through `Simulation::GetAgent` and down-casting.
  class World : public Object<World> {
  public:
    World() = default;

    /// @brief Build a `World` from a `world:` config block.
    ///
    /// Recognized keys: `frame` (reference frame, default `MOON_CI`); an optional
    /// `force_model:` sub-block (same schema as an `NBodyDynamics` `dynamics:`
    /// block: `bodies`, `integrator`, `dt`, `abstol`, `reltol`, `max_iter`,
    /// `autodiff`, and default `CR`/`area`/`mass`); an optional `gravity:` block
    /// (`body: MOON|EARTH|...` selecting the point-mass GM for `Gravity()`,
    /// default MOON); and an optional `dem:` block (`site_lat_deg`, `site_lon_deg`,
    /// `half_width_m`, `max_res_m`) loading a LOLA DEM and its ENU tangent frame.
    /// The epoch is taken from the global LuPNT epoch (`GetLupntEpoch()`), which
    /// `Simulation` sets from the scenario `epoch:`.
    explicit World(Config& world_config);

    /// @brief Register the owning `Simulation` (called by `Simulation::SetWorld`),
    /// enabling the `GetStateAt` truth facade to resolve agents by name.
    void SetSimulation(Simulation* sim) { sim_ = sim; }

    /// @brief TDB epoch [s past J2000] of simulation time t = 0.
    Real GetEpoch() const { return epoch_; }
    /// @brief Reference frame shared by the environment.
    Frame GetFrame() const { return frame_; }
    /// @brief The raw force-model config block, for inspection/serialization.
    const Config& GetForceModel() const { return force_model_; }
    /// @brief Whether a `force_model:` block was provided (so `MakeDynamics` works).
    bool HasForceModel() const { return has_force_model_; }

    /// @brief Construct a dynamics model for the shared force model, using the
    /// force model's default SRP spacecraft parameters (if any).
    Ptr<NBodyDynamics> MakeDynamics() const;

    /// @brief Construct a dynamics model for the shared force model, overriding
    /// the SRP spacecraft parameters (radiation-pressure coefficient, area
    /// [m^2], mass [kg]) for a specific spacecraft.
    Ptr<NBodyDynamics> MakeDynamics(double cr, double area_m2, double mass_kg) const;

    /// @brief Construct a **truth** dynamics model for the shared force model with
    /// the state-transition matrix disabled (`autodiff = false`) — truth
    /// propagation needs no STM. This is what a physical agent uses when it
    /// inherits the common `world:` force model instead of declaring its own
    /// `dynamics:` block. (The STM-carrying `MakeDynamics()` is for estimators.)
    Ptr<NBodyDynamics> MakeTruthDynamics() const;

    /// @brief The raw `plasma:` config block (ionosphere/plasmasphere signal-delay
    /// environment), if provided. This is a shared truth-environment property, so it
    /// lives under `world:` rather than in an individual receiver's application block;
    /// the GNSS ODTS app reads it from here. Requires `HasPlasma()`.
    const Config& GetPlasma() const { return plasma_; }
    /// @brief Whether a `plasma:` block was provided under `world:`.
    bool HasPlasma() const { return has_plasma_; }

    /// @brief Central-body gravitational parameter GM [m^3/s^2].
    double GetGM() const { return gm_; }
    /// @brief Point-mass central-body gravitational acceleration at position `r`
    /// (in the world frame) [m/s^2] — used by surface INS mechanization.
    Vec3d Gravity(const Vec3d& r) const;

    /// @brief Whether a `dem:` terrain block was provided.
    bool HasTerrain() const { return has_dem_; }
    /// @brief The loaded terrain DEM (for plotting/geometry). Requires `HasTerrain()`.
    const LunarDem& GetDem() const { return dem_; }
    /// @brief Terrain elevation at a local ENU offset (`east_m`, `north_m`) from the
    /// DEM/site center [m]. Requires `HasTerrain()`.
    double GetElevation(double east_m, double north_m) const;
    /// @brief Convert a local East-North-Up offset (from the site center) to a
    /// position in the world frame [m]. Requires `HasTerrain()`.
    Vec3d EnuToWorld(double east_m, double north_m, double up_m) const;
    /// @brief Rotation from the local ENU tangent frame to the world frame.
    const Mat3d& REnuToWorld() const { return R_enu2world_; }
    /// @brief Site center position in the world frame [m] (ENU origin).
    const Vec3d& SiteCenterWorld() const { return r_center_world_; }
    /// @brief Site latitude / east longitude [deg].
    double SiteLatDeg() const { return site_lat_deg_; }
    double SiteLonDeg() const { return site_lon_deg_; }

    /// @brief Read-only truth facade: the Cartesian state `[r; v]` of agent
    /// `agent_name` at simulation time `t` [s], in the agent's frame. Resolves the
    /// agent through the owning `Simulation` and calls its `GetStateAt` (which
    /// self-propagates without mutating the agent). Aborts if no `Simulation` is
    /// registered or the agent is unknown.
    Cart6 GetStateAt(const std::string& agent_name, Real t) const;

  private:
    Simulation* sim_ = nullptr;
    Real epoch_ = 0.0;
    Frame frame_ = Frame::MOON_CI;

    // Orbital force model.
    bool has_force_model_ = false;
    Config force_model_;  // NBodyDynamics-style force-model block

    // Ionosphere/plasmasphere signal-delay environment (shared truth property).
    bool has_plasma_ = false;
    Config plasma_;  // `plasma:` block; schema consumed by the GNSS ODTS app

    // Point-mass gravity (surface INS).
    double gm_ = 0.0;  // set from `gravity.body` (default MOON) in the constructor

    // Terrain / DEM (surface scenarios).
    bool has_dem_ = false;
    LunarDem dem_;
    double dem_cx_ = 0.0, dem_cy_ = 0.0;     // DEM native center coords
    Vec3d r_center_world_ = Vec3d::Zero();   // ENU origin in the world frame
    Mat3d R_enu2world_ = Mat3d::Identity();  // ENU -> world rotation
    double site_lat_deg_ = 0.0, site_lon_deg_ = 0.0;
  };

}  // namespace lupnt
