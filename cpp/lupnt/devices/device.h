#pragma once
#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>

#include "lupnt/core/config.h"
#include "lupnt/core/data_logger.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  class Agent;
  class Channel;

  /// @brief Base class for all onboard hardware devices (clocks, IMUs,
  /// cameras, comms radios) that can be attached to an `Agent`.
  ///
  /// Concrete devices (`Clock`, `Imu`, `Camera`, `Transmitter`/`Receiver`/
  /// `Transponder`, `GnssTransmitter`/`GnssReceiver`) derive from `Device`
  /// and override `Step`/`Setup`/`Log` to plug into the agent's simulation
  /// loop: `Setup` is called once during agent initialization to register the
  /// device (e.g. schedule periodic `Step` calls or bind to a `Channel`), and
  /// `Step` is then invoked by the `Simulation` event scheduler at the
  /// device's configured `frequency_` to advance its internal state/dynamics
  /// and produce measurements. `DataLogger::Log` is used by `Step` overrides
  /// to record device state to the simulation output.
  class Device : public DataLogger {
  public:
    Device() = default;

    /// @brief Construct a device from a YAML config node, reading the common
    /// `name` and (optional) `frequency` [Hz] fields used by `Agent`/
    /// `Simulation` to schedule periodic `Step` calls.
    Device(Config& config);
    virtual ~Device() = default;

    /// @brief Advance the device's internal state/dynamics by one simulation
    /// step at time `t` [s].
    ///
    /// Called by the `Simulation` event scheduler (via the periodic event
    /// registered in `Setup`) once per device update at the device's
    /// configured `frequency_`. Derived devices (`Clock`, `Imu`, `Camera`,
    /// comms radios) override this to propagate dynamics, generate
    /// measurements, and/or log state.
    ///
    /// @param t  Simulation time of this step [s]
    virtual void Step(Real t);

    /// @brief One-time initialization performed when the device is attached
    /// to an `Agent` and the simulation is built.
    ///
    /// The base implementation schedules a recurring `Step` event (priority
    /// `Event::Priority::DEVICE`) on the owning agent's `Simulation` if
    /// `frequency_ > 0`. Derived devices (e.g. `Transmitter`/`Receiver`/
    /// `Transponder`) override this to additionally resolve their `Channel`
    /// from the simulation by name.
    virtual void Setup();

    /// @brief Human-readable identifier for this device, used in log
    /// messages and `DataLogger` output paths (e.g. `"<name>/state/..."`).
    std::string GetName() const { return name_; }
    /// @brief Set the device's identifier (see `GetName`).
    void SetName(std::string name) { name_ = name; }

    /// @brief The `Agent` this device is attached to (e.g. a satellite,
    /// rover, or ground station), or `nullptr` if not yet attached.
    Agent* GetAgent() const { return agent_; }
    /// @brief Attach this device to an owning `Agent` (called by
    /// `Agent::AddDevice`).
    void SetAgent(Agent* agent) { agent_ = agent; }

    /// @brief Update rate at which `Step` is scheduled [Hz]; 0 means no
    /// periodic stepping is scheduled by the base `Setup`.
    Real GetFrequency() const { return frequency_; }
    /// @brief Set the periodic `Step` update rate [Hz] (see `GetFrequency`).
    void SetFrequency(Real frequency) { frequency_ = frequency; }

    /// @brief Record this device's current state to the simulation's
    /// `DataLogger` output at time `t` [s]. Overridden by derived devices
    /// (e.g. `Clock::Log`) to log model-specific quantities; the base
    /// implementation is a no-op.
    virtual void Log(Real time) { (void)time; };

  protected:
    Config config_;
    Real frequency_ = 0.0;
    std::string name_;
    Agent* agent_ = nullptr;
    Real time_ = 0.0;
  };

}  // namespace lupnt
