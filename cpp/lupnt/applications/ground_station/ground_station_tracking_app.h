#pragma once

#include <random>
#include <string>
#include <vector>

#include "lupnt/applications/application.h"
#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/config.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  class AgentWithDynamics;
  class GroundStationManagerApp;

  /// @brief Ground-station tracking (sensor) application.
  ///
  /// Runs on a `GroundStation` agent and generates that station's own two-way
  /// range and range-rate (Doppler) observations of a target satellite agent.
  /// On each scheduled `Step` it queries the target truth state (through the
  /// shared `World`), evaluates topocentric elevation, and -- if the target is
  /// above the elevation mask -- forms a noisy range/range-rate measurement and
  /// pushes it (with the station's own inertial state) to the
  /// `GroundStationManagerApp` on the manager agent, which runs the centralized
  /// orbit-determination filter over all stations.
  ///
  /// The station owns its observations (elevation and measurement series, for
  /// per-station visibility analysis); the manager owns the estimate.
  class GroundStationTrackingApp : public Application {
  public:
    GroundStationTrackingApp() = default;

    /// @brief Construct from the `application:` block of a `GroundStation` agent.
    /// Keys (besides base `name`/`frequency`): `target` (target agent name),
    /// `manager` (manager agent name), `elevation_mask_deg`, `use_range`,
    /// `use_range_rate`, `range_sigma_m`, `range_rate_sigma_mps`, `seed`.
    ///
    /// Optional signal-path / station-location corrections (all default off, so the
    /// observable stays a pure geometric range unless requested): `apply_solid_earth_tide`,
    /// `apply_troposphere` (+ `troposphere_pressure_hpa`, `troposphere_temperature_k`,
    /// `troposphere_humidity_pct`; unset pressure/temperature default to the standard
    /// atmosphere at the station height), `apply_ionosphere` (+ `ionosphere_vtec_tecu`,
    /// `signal_frequency_hz`), `apply_shapiro`. When enabled they are added to the *truth*
    /// observable this sensor generates; an estimator that does not model them then sees
    /// them as realistic errors. Path delays bias the range only; the solid Earth tide
    /// displaces the station position and so affects both range and range-rate.
    explicit GroundStationTrackingApp(Config& config);

    /// @brief Resolve the target agent, this station's geometry, and the manager
    /// app; register with the manager; schedule periodic measurement `Step`s.
    void Setup() override;

    /// @brief One measurement epoch: evaluate elevation, and (if visible) form and
    /// report a noisy range/range-rate observation to the manager.
    void Step(Real t) override;

    void Log(Real t) override;

    // Per-station observation series (valid after the run).
    const std::string& StationName() const { return station_name_; }
    const std::vector<double>& ElevationTime() const { return elev_t_; }         // [*] sim-rel [s]
    const std::vector<double>& ElevationDeg() const { return elev_deg_; }        // [*] deg
    const std::vector<double>& MeasurementTime() const { return meas_t_; }       // [*] sim-rel [s]
    const std::vector<double>& MeasurementRange() const { return meas_range_; }  // [*] m
    const std::vector<double>& MeasurementRangeRate() const { return meas_range_rate_; }  // [*] m/s

  protected:
    // Config
    std::string target_name_;
    std::string manager_name_;
    double elevation_mask_deg_ = 10.0;
    bool use_range_ = true;
    bool use_range_rate_ = true;
    double range_sigma_m_ = 10.0;
    double range_rate_sigma_mps_ = 1.0e-3;
    int seed_ = 42;

    // Optional truth-model corrections (all default off).
    bool apply_solid_earth_tide_ = false;
    bool apply_troposphere_ = false;
    bool apply_ionosphere_ = false;
    bool apply_shapiro_ = false;
    double tropo_pressure_hpa_ = -1.0;    // <=0 -> standard atmosphere at station height
    double tropo_temperature_k_ = -1.0;   // <=0 -> standard atmosphere at station height
    double tropo_humidity_pct_ = 50.0;    // relative humidity [%]
    double iono_vtec_tecu_ = 10.0;        // vertical TEC [TECU]
    double signal_frequency_hz_ = 8.4e9;  // carrier for the ionosphere delay (X-band)

    // Resolved at Setup
    GroundStationManagerApp* manager_ = nullptr;
    int station_id_ = 0;
    std::string station_name_;
    Real epoch0_ = 0.0;              // absolute TDB epoch of sim time t = 0
    Vec3 station_r_ = Vec3::Zero();  // station position, station body-fixed frame
    Frame station_frame_ = Frame::UNDEFINED;
    double station_lat_rad_ = 0.0;      // geodetic latitude (for tropo/tide) [rad]
    double station_height_m_ = 0.0;     // height above the ellipsoid [m]
    bool corrections_enabled_ = false;  // any of the optional corrections requested
    std::mt19937 noise_rng_;

    // Observation series
    std::vector<double> elev_t_, elev_deg_;
    std::vector<double> meas_t_, meas_range_, meas_range_rate_;
  };

}  // namespace lupnt
