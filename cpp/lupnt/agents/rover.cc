#include "lupnt/agents/rover.h"

#include "lupnt/core/logger.h"

namespace lupnt {

  Rover::Rover(Config& config) : AgentWithDynamics() {
    config_ = config;
    if (config["name"])
      SetName(config["name"].as<std::string>());
    else
      SetName(GetId());
    if (config["frequency"]) SetFrequency(config["frequency"].as<Real>());
    Logger::Debug(fmt::format("Creating Rover {}", GetName()), "Rover");

    // Thin agent: no dynamics block required; the hosted application drives the truth state.
    // Attach the onboard navigation application (e.g. SurfaceRoverNavApp) from the config.
    CreateApplication(config);
  }

  REGISTER_FACTORY_CLASS(Agent, Rover)

}  // namespace lupnt
