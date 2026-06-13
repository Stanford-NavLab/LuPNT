#pragma once

#include "lupnt/core/config.h"
#include "lupnt/devices/device.h"

namespace lupnt {

  /// @brief Onboard imaging camera device: holds the camera's pose
  /// (position/orientation in the host body frame) and image dimensions for
  /// surface-imaging operations.
  ///
  /// Attached to a rover/lander/orbiter `Agent` to model an imaging payload
  /// used by surface-operations applications (e.g. crater/feature imaging in
  /// `applications/`); rendering/measurement generation from the camera pose
  /// is not yet implemented (see commented-out `renderer_` member).
  class Camera : public Device {
  public:
    /// @brief Construct a camera device from a YAML config node (see
    /// `Device::Device(Config&)`).
    Camera(Config& config);
    /// @brief Camera device step; currently a no-op placeholder. Overrides
    /// `Device::Step`.
    void Step(Real t) override;
    /// @brief Camera device setup; currently a no-op placeholder. Overrides
    /// `Device::Setup`.
    void Setup() override;

  protected:
    Vec3 position_ = Vec3::Zero();
    Mat3 orientation_ = Mat3::Identity();
    int width_ = 0;
    int height_ = 0;

    // Renderer* renderer_ = nullptr;
  };

}  // namespace lupnt
