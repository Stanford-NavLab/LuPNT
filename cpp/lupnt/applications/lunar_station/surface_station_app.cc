#include "lupnt/applications/lunar_station/surface_station_app.h"

#include "lupnt/core/asset_factory.h"
#include "lupnt/core/logger.h"

namespace lupnt {

  // SurfaceStationApp
  SurfaceStationApp::SurfaceStationApp(Config& config) : Application(config) {
    Logger::Debug(fmt::format("Creating {}", name_), "SurfaceStationApp");
  }
  void SurfaceStationApp::Setup() {
    Logger::Debug(fmt::format("Setting up {}", name_), "SurfaceStationApp");
  }
  void SurfaceStationApp::Step(Real t) { Logger::Debug("Step", "SurfaceStationApp", t); }

  REGISTER_FACTORY_CLASS(Application, SurfaceStationApp)

}  // namespace lupnt
