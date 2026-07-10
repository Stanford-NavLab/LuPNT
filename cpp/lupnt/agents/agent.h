/**
 * @file agent.h
 * @author Stanford NAV LAB
 * @brief List of agents
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <memory>

#include "lupnt/core/asset_factory.h"
#include "lupnt/core/config.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/data_logger.h"
#include "lupnt/core/object.h"
#include "lupnt/dynamics/attitude_dynamics.h"
#include "lupnt/dynamics/dynamics.h"
#include "lupnt/states/state.h"

namespace lupnt {

  class Simulation;
  class World;
  class Channel;
  class Application;
  class Device;
  class Agent;

  using DeviceFactory = AssetFactory<Device, YAML::Node&>;
  using AgentFactory = AssetFactory<Agent, YAML::Node&>;
  using ApplicationFactory = AssetFactory<Application, YAML::Node&>;
  using DynamicsFactory = AssetFactory<Dynamics, YAML::Node&>;
  using ChannelFactory = AssetFactory<Channel, YAML::Node&>;

  class Agent : public Object<Agent>, public DataLogger {
  protected:
    // Config
    std::string name_;
    Simulation* sim_ = nullptr;
    Config config_;
    Real frequency_ = 0.0;
    Real time_ = 0.0;
    State state_ = State();

    // Devices
    std::map<std::string, Ptr<Device>> devices_;

    // Application
    Ptr<Application> application_;

  public:
    Agent() = default;

    /// @brief Construct an agent from a YAML config node, creating and
    /// attaching its devices (and, if specified, application).
    ///
    /// Called by `AssetFactory<Agent, Config&>` (via `AgentFactory`) when the
    /// simulation builder instantiates the `agents:` section of a scenario
    /// config; reads `name`, `frequency`, `devices`, and `application` keys
    /// and creates each listed device/application through its own factory,
    /// wiring it back to this agent via `AddDevice`/`SetApplication`.
    /// @param agent_config YAML config node for this agent
    Agent(Config& agent_config);

    /// @brief Get the `Simulation` this agent is registered with.
    /// @return Pointer to the owning `Simulation` (aborts via `LUPNT_CHECK`
    ///         if the agent has not yet been added to one).
    Simulation* GetSimulation() const {
      LUPNT_CHECK(sim_, "Simulation not set", "Agent");
      return sim_;
    }
    /// @brief Register this agent with a `Simulation`, called by
    /// `Simulation::AddAgent` so that `GetSimulation`/`Setup` can access the
    /// event scheduler and Cesium viewer.
    void SetSimulation(Simulation* sim) { sim_ = sim; }

    /// @brief Get the shared `World` (read-only environment) from the owning
    /// simulation, or nullptr if the scenario defined no `world:` block. Returns
    /// the same object as `GetSimulation()->GetWorld()`.
    World* GetWorld() const;

    /// @brief Get the agent's current simulation time [s].
    Real GetTime() const { return time_; }
    /// @brief Set the agent's current simulation time [s].
    void SetTime(Real t) { time_ = t; }

    /// @brief Get the agent's name (used as a prefix for device names and
    /// data-logger keys).
    std::string GetName() const { return name_; }
    /// @brief Set the agent's name.
    void SetName(std::string name) { name_ = name; }

    /// @brief Get the map of devices (camera, clock, radio, etc.) attached to
    /// this agent, keyed by device name.
    std::map<std::string, Ptr<Device>> GetDevices() const { return devices_; }

    /// @brief Attach a `Device` to this agent, called during construction
    /// (from the `devices:` section of the config) or by application setup
    /// code that adds devices programmatically. Sets the device's back-link
    /// to this agent and registers it under `device->GetName()`.
    /// @param device Device to attach; ownership is transferred into
    ///        `devices_`
    void AddDevice(Ptr<Device> device);

    /// @brief Look up a previously-attached device by name, trying both the
    /// bare name and `<agent_name>/<name>` (the convention used by
    /// `AddDevice`/`Agent::Agent`). Used by application/measurement code to
    /// retrieve e.g. a `Clock` or `Transmitter` device from its owning agent.
    /// Aborts via `LUPNT_CHECK` if no matching device is found.
    /// @param name Device name (with or without the `<agent_name>/` prefix)
    /// @return     The matching device
    Ptr<Device> GetDevice(const std::string& name) const;

    /// @brief Create and attach this agent's `Application` from the
    /// `application:` block of a config node (if present). Factored out of the
    /// `Agent(Config&)` constructor so agents that don't chain to it (e.g.
    /// `GroundStation`, which builds itself via `ConfigureFromConfig`) can still
    /// opt into an application. Requires `name_` to already be set.
    void CreateApplication(Config& config);

    /// @brief Get the `Application` attached to this agent (e.g. an LNSS,
    /// rover, or surface-station application), if any.
    Ptr<Application> GetApplication() const { return application_; }
    /// @brief Attach an `Application` to this agent, called during
    /// construction (from the `application:` section of the config). Wires the
    /// application's back-pointer to this agent so it can reach the owning
    /// agent/simulation from its `Setup`/`Step`/`Log`.
    void SetApplication(Ptr<Application> app);

    /// @brief Perform one-time setup before the simulation starts: schedules
    /// a periodic `Step` callback with the simulation event scheduler (if
    /// `frequency_ > 0`) and calls `Setup` on each attached device.
    ///
    /// Called once by `Simulation::Setup` for every registered agent, before
    /// the main time-stepping loop begins. `AgentWithDynamics` and
    /// `GnssConstellation` override/extend this for their own initialization
    /// needs.
    virtual void Setup();

    /// @brief Advance the agent by one simulation step at time `time`
    /// (propagate dynamics, update devices/application, and log state).
    ///
    /// Invoked periodically by the simulation's event scheduler at the rate
    /// configured via `frequency_` (see `Setup`). The base implementation is
    /// a no-op (logging only); `AgentWithDynamics::Step` overrides it to call
    /// `Propagate` and `Log`.
    /// @param time Current simulation time [s]
    virtual void Step(Real time);

    /// @brief Get the agent's current state. The base `Agent` class has no
    /// intrinsic state representation and always returns an empty `State()`;
    /// derived classes such as `AgentWithDynamics` override this to return
    /// their Cartesian/orbital state.
    State GetState() const { return State(); }
    /// @brief Set the agent's current state (base-class no-op storage; see
    /// `AgentWithDynamics::SetState` for the overload that actually drives
    /// dynamics propagation).
    void SetState(const State& state) { state_ = state; }

    /// @brief Compute (or interpolate/propagate) the agent's Cartesian
    /// state `[r; v]` at an arbitrary time `t`, without mutating the agent's
    /// stored state/time.
    ///
    /// Pure virtual: every concrete agent (Satellite, Rover, GroundStation,
    /// SurfaceStation, ...) must provide this so that measurement models and
    /// other agents can query "where was/will this agent be at time `t`"
    /// (e.g. for light-time-corrected GNSS pseudorange computation).
    /// @param t Query time [s]
    /// @return  Cartesian state `[r; v]` at time `t`, in the agent's
    ///          dynamics frame
    virtual Cart6 GetStateAt(Real t) const = 0;

    /// @brief Write this agent's state/devices/application to the
    /// `DataLogger`. Called once per `Step` (after propagation) by
    /// `AgentWithDynamics::Step`; the base implementation is a no-op.
    /// @param time Current simulation time [s]
    virtual void Log(Real time);

    /// @brief Push this agent's current state/geometry to the Cesium 3D
    /// viewer (if one is attached to the simulation). Called once per `Step`
    /// alongside `Log`; the base implementation is a no-op.
    virtual void LogCesium();
  };

  class AgentWithDynamics : public Agent {
  protected:
    // Config
    BodyId body_id_;
    Real precompute_ = 0.0;  // [s]
    Real frequency_ = 0.0;   // [Hz]

    // State
    Real time_ = 0.0;  // [s]
    Cart6 state_;
    Attitude attitude_;
    State control_;

    // Dynamics
    Ptr<Dynamics> dynamics_;
    Ptr<AttitudeDynamics> attitude_dynamics_;

  public:
    AgentWithDynamics() = default;

    /// @brief Construct an agent with orbital and attitude dynamics from a
    /// YAML config node.
    ///
    /// Extends `Agent::Agent`: in addition to devices/application, creates
    /// the `dynamics:` and `attitude_dynamics:` blocks via their respective
    /// `AssetFactory` (defaulting attitude dynamics to
    /// `FixedAttitudeDynamics` if not specified), initializes `state_`/
    /// `attitude_` to zero in the `MOON_CI` frame, and reads the optional
    /// `precompute` interval. Used as the base constructor for `Satellite`,
    /// `Rover`, `GroundStation`, and `SurfaceStation`.
    /// @param agent_config YAML config node for this agent
    AgentWithDynamics(Config& agent_config);

    /// @brief Base-class setup (schedule periodic `Step`, set up devices);
    /// `AgentWithDynamics` currently adds no further initialization beyond
    /// `Agent::Setup`.
    virtual void Setup() override;

    /// @brief Advance this agent by one simulation step: propagate
    /// orbital/attitude dynamics to `time` via `Propagate`, then log the
    /// resulting state via `Log`.
    /// @param time Target simulation time [s]
    virtual void Step(Real time) override;

    /// @brief Get the agent's current simulation time [s].
    Real GetTime() const { return time_; }
    /// @brief Set the agent's current simulation time [s] (without
    /// propagating; use `Propagate` to advance the state consistently).
    void SetTime(Real t) { time_ = t; }

    /// @brief Get the agent's `Step` scheduling frequency [Hz].
    Real GetFrequency() const { return frequency_; }
    /// @brief Set the agent's `Step` scheduling frequency [Hz], used by
    /// `Setup` to register the periodic simulation callback.
    void SetFrequency(Real frequency) { frequency_ = frequency; }

    /// @brief Get the agent's current Cartesian state `[r; v]` (in the frame
    /// returned by `dynamics_->GetFrame()`, typically `MOON_CI`).
    State GetState() const { return state_; }
    /// @brief Set the agent's current Cartesian state `[r; v]`, e.g. to
    /// initialize a scenario or apply a filter/estimator correction.
    void SetState(const State& x) { state_ = x; }

    /// @brief Get the agent's current attitude (quaternion + angular rate).
    Attitude GetAttitude() const { return attitude_; }
    /// @brief Set the agent's current attitude (quaternion + angular rate).
    void SetAttitude(const Attitude& attitude) { attitude_ = attitude; }

    /// @brief Get the agent's current control input (e.g. rover
    /// linear/angular velocity commands), passed to `dynamics_->Propagate`.
    State GetControl() const { return control_; }
    /// @brief Set the agent's current control input, applied on the next
    /// call to `Propagate`/`Step`.
    void SetControl(const State& control) { control_ = control; }

    /// @brief Get the orbital/translational `Dynamics` model driving
    /// `Propagate`.
    /// @return Raw pointer to the dynamics model (aborts via `LUPNT_CHECK`
    ///         if none has been set).
    Dynamics* GetDynamics() const {
      LUPNT_CHECK(dynamics_, "Dynamics not set", "Agent");
      return dynamics_.get();
    }
    /// @brief Set the orbital/translational dynamics model, called during
    /// construction (from the `dynamics:` config block) or to swap dynamics
    /// models at runtime.
    void SetDynamics(Ptr<Dynamics> dyn) { dynamics_ = std::move(dyn); }

    /// @brief Get the `AttitudeDynamics` model driving attitude propagation
    /// in `Propagate`/`GetAttitudeAt`.
    /// @return Raw pointer to the attitude dynamics model (aborts via
    ///         `LUPNT_CHECK` if none has been set).
    AttitudeDynamics* GetAttitudeDynamics() const {
      LUPNT_CHECK(attitude_dynamics_, "Attitude dynamics not set", "Agent");
      return attitude_dynamics_.get();
    }
    /// @brief Set the attitude dynamics model, called during construction
    /// (from the `attitude_dynamics:` config block, defaulting to
    /// `FixedAttitudeDynamics`) or to swap models at runtime.
    void SetAttitudeDynamics(Ptr<AttitudeDynamics> dyn) { attitude_dynamics_ = std::move(dyn); }

    /// @brief Set the ID of the central body this agent's dynamics are
    /// referenced to (e.g. `BodyId::EARTH`, `BodyId::MOON`).
    void SetBodyId(BodyId body_id) { body_id_ = body_id; }
    /// @brief Get the ID of the central body this agent's dynamics are
    /// referenced to.
    BodyId GetBodyId() const { return body_id_; }

    /// @brief Advance the agent's stored orbital state, attitude, and time
    /// from `time_` to `t` by integrating `dynamics_`/`attitude_dynamics_`.
    ///
    /// Called once per `Step` by `AgentWithDynamics::Step` (and overridden by
    /// `Rover::Propagate` for its 2D surface-dynamics model). Mutates
    /// `state_`, `attitude_`, and `time_` in place using the current
    /// `control_`.
    /// @param t Target simulation time [s]
    virtual void Propagate(Real t);

    /// @brief Compute the agent's Cartesian state `[r; v]` at time `t`
    /// without mutating the stored state.
    ///
    /// If `t` matches the agent's current `time_` (within `EPS`), returns
    /// the cached `state_` directly; otherwise propagates a copy of `state_`
    /// from `time_` to `t` via `dynamics_->Propagate`. Used by measurement
    /// models and other agents to query this agent's position/velocity at an
    /// arbitrary (e.g. light-time-corrected) epoch.
    /// @param t Query time [s]
    /// @return  Cartesian state `[r; v]` at time `t`
    virtual Cart6 GetStateAt(Real t) const override;

    /// @brief Compute the agent's attitude at time `t` without mutating the
    /// stored attitude.
    ///
    /// If `t` matches `time_` (within `EPS`), or no attitude dynamics are
    /// set, returns the cached `attitude_` directly; otherwise propagates a
    /// copy of `attitude_` from `time_` to `t` via
    /// `attitude_dynamics_->Propagate`.
    /// @param t Query time [s]
    /// @return  Attitude (quaternion + angular rate) at time `t`
    virtual Attitude GetAttitudeAt(Real t) const;

    // Log

    /// @brief Log this agent's state, attitude, and control vectors (and the
    /// attached application's state, if any) to the `DataLogger`. Called
    /// once per `Step` by `AgentWithDynamics::Step`.
    /// @param time Current simulation time [s]
    virtual void Log(Real time) override;

    /// @brief Push this agent's geometry to the Cesium 3D viewer. Currently a
    /// no-op for the base `AgentWithDynamics`; `Satellite`/`SurfaceStation`
    /// override it (see their respective headers).
    virtual void LogCesium() override;
  };

};  // namespace lupnt
