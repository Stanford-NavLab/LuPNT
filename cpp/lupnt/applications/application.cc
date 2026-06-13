#include "lupnt/applications/application.h"

#include "lupnt/agents/agent.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/core/error.h"
#include "lupnt/core/logger.h"
#include "lupnt/core/simulation.h"
#include "lupnt/devices/device.h"
#include "lupnt/interfaces/yaml.h"

namespace lupnt {

  // Application
  Application::Application(Config& config) {
    name_ = config["name"] ? config["name"].as<std::string>() : GetId();
    config_ = config;
    Logger::Debug(fmt::format("Creating {}", name_), "Application");

    if (config["frequency"]) frequency_ = config["frequency"].as<Real>();
  }

  void Application::Setup() {
    Logger::Debug(fmt::format("Setting up {}", name_), "Application");

    // Schedule
    if (frequency_ > 0.0) {
      agent_->GetSimulation()->Schedule(
          0.0, [this](Real t) { Step(t); }, frequency_, Event::Priority::APPLICATION);
      Logger::Info(fmt::format("Scheduled {}, t={} s, f={} Hz", name_, 0.0, frequency_), "Agent");
    } else {
      Logger::Warn(fmt::format("No frequency set for {}", name_), "Agent");
    }
  }

  void Application::Log(Real t) {
    Logger::Debug(fmt::format("Logging {}", name_), "Application", t);
  }

  // Define the GetRegistry function for this specialization (must come before explicit
  // instantiation)
  template <> std::unordered_map<std::string, AssetFactory<Application, Config&>::Creator>&
  AssetFactory<Application, Config&>::GetRegistry() {
    return Registry();
  }

  // Explicit template instantiation to ensure single registry across library boundaries
  template class AssetFactory<Application, Config&>;

}  // namespace lupnt
