#include "lupnt/agents/surface_station.h"

#include "lupnt/conversions/frame_conversions.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/core/logger.h"
#include "lupnt/interfaces/yaml.h"

namespace lupnt {

  // SurfaceStation
  SurfaceStation::SurfaceStation(Config& config) : AgentWithDynamics(config) {
    Logger::Debug(fmt::format("Creating SurfaceStation {}", name_), "SurfaceStation");
  }

  void SurfaceStation::Log(Real time) { Agent::Log(time); }

  void SurfaceStation::LogCesium() {}

  REGISTER_FACTORY_CLASS(Agent, SurfaceStation)

}  // namespace lupnt
