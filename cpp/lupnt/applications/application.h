#pragma once

#include <yaml-cpp/yaml.h>

#include "lupnt/core/config.h"
#include "lupnt/core/definitions.h"
#include "lupnt/core/object.h"

namespace lupnt {

  class Agent;

  /// @brief Base class for top-level simulation "application" logic attached to an Agent
  /// (e.g. LanderNavApp, SurfaceStationApp) -- the part of an agent that runs its navigation
  /// filter / mission logic on a periodic schedule.
  ///
  /// An Application is created and owned by an Agent (see the Agent(Config&) constructor,
  /// which builds it from the `application:` block of the agent's Config and calls
  /// `SetApplication(...)`). Derived classes wire together the agent's dynamics,
  /// devices, and a Filter (EKF/UKF/etc.) into a runnable scenario: Setup() initializes the
  /// filter state/covariance and schedules periodic Step() calls, Step() performs one
  /// predict/update cycle, and Log() records results each step (Agent::Log forwards to it).
  class Application : public Object<Application> {
  public:
    Application() = default;

    /// @brief Construct an Application from a YAML Config node.
    ///
    /// Called via AssetFactory<Application, Config&>::Create from the Agent(Config&)
    /// constructor when an agent's config contains an `application:` block. Reads the
    /// optional `name` (defaults to GetId()) and `frequency` [Hz] entries and stores the
    /// full config for use by derived-class constructors (e.g. LanderNavApp reads
    /// `dynamics`/`filter` sub-configs from it).
    ///
    /// @param config YAML configuration node for this application (must outlive derived
    ///               objects that retain a reference into it).
    Application(Config& config);
    virtual ~Application() = default;

    /// @brief Initialize the application and schedule its periodic Step() calls.
    ///
    /// Base implementation: if GetFrequency() > 0, schedules `Step(t)` on the owning
    /// agent's Simulation at Event::Priority::APPLICATION, starting at t=0 with period
    /// 1/GetFrequency() seconds; otherwise logs a warning that no frequency is set.
    /// Derived classes (e.g. LanderNavApp::Setup, SurfaceStationApp::Setup) override this to
    /// additionally initialize their Filter's time/state/covariance and dynamics/process
    /// noise callbacks before (optionally) calling the base behavior.
    virtual void Setup();

    /// @brief Run one simulation step of this application's mission logic at time `t`.
    ///
    /// Pure virtual: derived classes implement the actual filter predict/update cycle
    /// here (e.g. LanderNavApp::Step computes a control input and logs state; the LNSS
    /// application's Step would run the GNSS measurement-update). Invoked periodically
    /// by the Simulation event scheduled in Setup(), at the application's configured
    /// frequency.
    ///
    /// @param t Current simulation time [s, since simulation epoch]
    virtual void Step(Real t) = 0;

    /// @brief Log this application's current state/diagnostics to the DataLogger.
    ///
    /// Base implementation only emits a debug message. Called once after Setup() and
    /// thereafter from Agent::Log(t) (which forwards to `application_->Log(time)` each
    /// time the owning agent logs). Derived classes (e.g. LanderNavApp::Log) override this to
    /// additionally log their Filter's state/covariance and any error metrics.
    ///
    /// @param t Current simulation time [s, since simulation epoch]
    virtual void Log(Real t);

    /// @brief Get the application's name (used as a prefix for logging/data keys).
    std::string GetName() const { return name_; }
    /// @brief Set the application's name.
    void SetName(std::string name) { name_ = name; }

    /// @brief Get the configured Step() call frequency [Hz].
    Real GetFrequency() const { return frequency_; }
    /// @brief Set the Step() call frequency [Hz]; used by Setup() to schedule periodic steps.
    void SetFrequency(Real frequency) { frequency_ = frequency; }

    /// @brief Get the Agent that owns this application.
    Agent* GetAgent() const { return agent_; }
    /// @brief Set the Agent that owns this application (called by Agent::SetApplication).
    void SetAgent(Agent* agent) { agent_ = agent; }

  protected:
    std::string name_;
    Real frequency_ = 0.0;
    Agent* agent_ = nullptr;
    Config config_;
  };

};  // namespace lupnt
