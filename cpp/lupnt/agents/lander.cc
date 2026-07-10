#include "lupnt/agents/lander.h"

#include "lupnt/core/logger.h"

namespace lupnt {

  Lander::Lander(Config& config) : AgentWithDynamics() {
    config_ = config;
    if (config["name"])
      SetName(config["name"].as<std::string>());
    else
      SetName(GetId());
    if (config["frequency"]) SetFrequency(config["frequency"].as<Real>());
    Logger::Debug(fmt::format("Creating Lander {}", GetName()), "Lander");

    // Thin agent: no dynamics block required; the hosted application drives the truth state.
    CreateApplication(config);
  }

  REGISTER_FACTORY_CLASS(Agent, Lander)

}  // namespace lupnt
