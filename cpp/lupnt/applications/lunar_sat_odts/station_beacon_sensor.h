#pragma once

#include <random>
#include <string>
#include <vector>

#include "lupnt/applications/application.h"
#include "lupnt/core/config.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  class Spacecraft;
  class GroundOdtsApp;

  /// @brief Ground-tracking *sensor*, hosted on a physical `SurfaceStation` agent.
  ///
  /// Each epoch it synthesizes the one-way pseudorange (and optional Doppler) to every visible
  /// satellite from the satellites' truth and its own fixed lunar-surface position, then pushes
  /// those observations to the `GroundOdtsApp` estimator on the `SurfaceStationManager`
  /// (`manager_->AddMeasurement(...)`). This is the ISL-scenario analog of the ground-station
  /// example's `GroundStationTrackingApp`: the sensor knows the truth (it is the hardware that
  /// observes the physical world); the estimator it feeds does not.
  class StationBeaconSensor : public Application {
  public:
    StationBeaconSensor() = default;
    explicit StationBeaconSensor(Config& config);

    void Setup() override;
    void Step(Real t) override;
    void Log(Real /*t*/) override {}

  protected:
    // Config
    std::string manager_name_;
    std::vector<std::string> sat_names_;
    double lat_deg_ = 0.0, lon_deg_ = 0.0, alt_m_ = 0.0;
    double elevation_mask_deg_ = 5.0;
    double pseudorange_sigma_m_ = 1.0;
    bool include_doppler_ = true;
    double doppler_sigma_mps_ = 1.0e-3;
    int seed_ = 42;

    // Resolved
    bool initialized_ = false;
    GroundOdtsApp* manager_ = nullptr;
    std::vector<Spacecraft*> sats_;
    std::vector<int> sat_index_;  // filter block index per satellite
    Vec3 r_bf_;                   // station position in MOON_PA [m]
    std::mt19937 rng_;

    void Initialize();
  };

}  // namespace lupnt
