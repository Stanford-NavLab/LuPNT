#include "lupnt/devices/camera.h"

#include <fmt/format.h>

#include "lupnt/core/asset_factory.h"
#include "lupnt/core/logger.h"
#include "lupnt/devices/device.h"
#include "lupnt/measurements/measurement.h"
#include "lupnt/measurements/measurement_utils.h"

namespace lupnt {

  // Constructor
  Camera::Camera(Config& config) : Device(config) {
    Logger::Debug(fmt::format("Creating {}", name_), "Camera");
  }

  // Setup
  void Camera::Setup() { Logger::Debug(fmt::format("Setting up {}", name_), "Camera"); }

  // Step
  void Camera::Step(Real t) { Logger::Debug("Step", "Camera", t); }

  REGISTER_FACTORY_CLASS(Device, Camera)

}  // namespace lupnt
