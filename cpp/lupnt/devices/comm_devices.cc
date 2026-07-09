#include "lupnt/devices/comm_devices.h"

#include <fmt/format.h>

#include "lupnt/core/asset_factory.h"
#include "lupnt/core/logger.h"
#include "lupnt/devices/device.h"
#include "lupnt/measurements/channel.h"
#include "lupnt/simulations/simulation.h"

namespace lupnt {

  // Transmitter **************************************************************

  // Constructor
  Transmitter::Transmitter(Config& config) : Device(config) {
    std::string device_name = config["name"] ? config["name"].as<std::string>() : "transmitter";
    Logger::Debug("Creating", "Transmitter");
  }

  // Setup
  void Transmitter::Setup() {
    Logger::Debug(fmt::format("Setting up {}", name_), "Transmitter");
    if (config_["channel"]) {
      auto sim = agent_->GetSimulation();
      channel_ = sim->GetChannel(config_["channel"].as<std::string>());
      Logger::Debug(fmt::format("Channel {} set for {}", channel_->GetName(), name_),
                    "Transmitter");
    } else {
      Logger::Warn(fmt::format("No channel set for {}", name_), "Transmitter");
    }
  }

  // Step
  void Transmitter::Step(Real t) { Logger::Debug("Step", "Transmitter", t); }

  // Called by channel requesting data
  void Transmitter::Send(Real t) {
    LUPNT_CHECK(channel_, "Channel not set", "Transmitter");
    Logger::Debug(fmt::format("Transmitter {} sending", name_), "Transmitter", t);
  }

  // Called by external device providing data
  void Transmitter::Send(Real t, void* data) {
    LUPNT_CHECK(channel_, "Channel not set", "Transmitter");
    LUPNT_CHECK(data, "Data not set", "Transmitter");
    Logger::Debug(fmt::format("Transmitter {} sending", name_), "Transmitter", t);
    channel_->Send(this, t, data);
  }

  REGISTER_FACTORY_CLASS(Device, Transmitter)

  // Receiver **************************************************************

  // Constructor
  Receiver::Receiver(Config& config) : Device(config) {
    Logger::Debug(fmt::format("Creating {}", name_), "Receiver");
  }
  void Receiver::Setup() {
    Logger::Debug(fmt::format("Setting up {}", name_), "Receiver");
    if (config_["channel"]) {
      LUPNT_CHECK(agent_, "Agent not set", "Receiver");
      auto sim = agent_->GetSimulation();
      LUPNT_CHECK(sim, "Simulation not set", "Receiver");
      channel_ = sim->GetChannel(config_["channel"].as<std::string>());
      Logger::Debug(fmt::format("Channel {} set for {}", channel_->GetName(), name_), "Receiver");
    } else {
      Logger::Warn(fmt::format("No channel set for {}", name_), "Receiver");
    }
  }

  void Receiver::Step(Real t) { Logger::Debug("Step", "Receiver", t); }

  // Called by channel requesting data
  void Receiver::Receive(Real t) {
    Logger::Debug(fmt::format("Receiver {} receiving", name_), "Receiver", t);
  }

  // Called by channel providing data
  void Receiver::Receive(Real t, void* data) {
    Logger::Debug(fmt::format("Receiver {} receiving", name_), "Receiver", t);
    data_.push_back(data);
  }

  REGISTER_FACTORY_CLASS(Device, Receiver)

  // Transponder **************************************************************

  // Constructor
  Transponder::Transponder(Config& config) : Device(config) {
    Logger::Debug(fmt::format("Creating {}", name_), "Transponder");
  }

  // Setup
  void Transponder::Setup() {
    Logger::Debug(fmt::format("Setting up {}", name_), "Transponder");
    if (config_["channel"]) {
      LUPNT_CHECK(agent_, "Agent not set", "Transponder");
      auto sim = agent_->GetSimulation();
      LUPNT_CHECK(sim, "Simulation not set", "Transponder");
      channel_ = sim->GetChannel(config_["channel"].as<std::string>());
      Logger::Debug(fmt::format("Channel {} set for {}", channel_->GetName(), name_),
                    "Transponder");
    } else {
      Logger::Warn(fmt::format("No channel set for {}", name_), "Transponder");
    }
  }

  // Step
  void Transponder::Step(Real t) { Logger::Debug("Step", "Transponder", t); }

  // Called by channel requesting data
  void Transponder::Send(Real t) {
    LUPNT_CHECK(channel_, "Channel not set", "Transponder");
    Logger::Debug(fmt::format("Transponder {} sending", name_), "Transponder", t);
  }

  // Called by external device providing data
  void Transponder::Send(Real t, void* data) {
    LUPNT_CHECK(channel_, "Channel not set", "Transponder");
    LUPNT_CHECK(data, "Data not set", "Transponder");
    Logger::Debug(fmt::format("Transponder {} sending", name_), "Transponder", t);
    channel_->Send(this, t, data);
  }

  // Called by external device requesting data
  void Transponder::Receive(Real t) {
    LUPNT_CHECK(channel_, "Channel not set", "Transponder");
    Logger::Debug(fmt::format("Transponder {} receiving", name_), "Transponder", t);
    channel_->Receive(this, t);
  }

  // Called by channel providing data
  void Transponder::Receive(Real t, void* data) {
    LUPNT_CHECK(data, "Data not set", "Transponder");
    Logger::Debug(fmt::format("Transponder {} receiving", name_), "Transponder", t);
    data_.push_back(data);
  }

  REGISTER_FACTORY_CLASS(Device, Transponder)

}  // namespace lupnt
