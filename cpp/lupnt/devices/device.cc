#include "lupnt/devices/device.h"

#include <fmt/format.h>

#include "lupnt/agents/agent.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/core/definitions.h"
#include "lupnt/core/event.h"
#include "lupnt/core/logger.h"
#include "lupnt/measurements/channel.h"
#include "lupnt/simulations/simulation.h"

namespace lupnt {

  Device::Device(Config& config) {
    config_ = config;
    name_ = config["name"] ? config["name"].as<std::string>() : "device";
    Logger::Debug(fmt::format("Creating {}", name_), "Device");
    if (config["frequency"]) {
      frequency_ = config["frequency"].as<double>();
    }
  }
  void Device::Step(Real t) { Logger::Debug(fmt::format("Step {}", name_), "Device", t); }

  void Device::Setup() {
    // Step
    if (frequency_ > 0.0) {
      agent_->GetSimulation()->Schedule(
          0.0, [this](Real t) { Step(t); }, frequency_, Event::Priority::DEVICE);
      Logger::Debug(fmt::format("Scheduled {}, t={} s, f={} Hz", name_, 0.0, frequency_), "Device");
    }
  }

  // Define the GetRegistry function for this specialization (must come before explicit
  // instantiation)
  template <> std::unordered_map<std::string, AssetFactory<Device, Config&>::Creator>&
  AssetFactory<Device, Config&>::GetRegistry() {
    return Registry();
  }

  // Explicit template instantiation to ensure single registry across library boundaries
  template class AssetFactory<Device, Config&>;

}  // namespace lupnt
