#pragma once

#include <map>
#include <string>
#include <vector>

#include "lupnt/applications/application.h"
#include "lupnt/core/config.h"
#include "lupnt/core/definitions.h"
#include "lupnt/numerics/filters/ekf.h"

namespace lupnt {

  class Spacecraft;

  /// @brief One station->satellite pseudorange (+ optional Doppler) observation, pushed to
  /// the `GroundOdtsApp` estimator by a `StationBeaconSensor` on a physical `SurfaceStation`.
  struct IslStationObs {
    double t = 0.0;
    int sat_index = 0;      // block index of the observed satellite in the filter state
    Vec3d station_mci;      // station position in MOON_CI [m]
    Vec3d station_vel_mci;  // station velocity in MOON_CI [m/s]
    double pseudorange_m = 0.0;
    bool has_doppler = false;
    double doppler_mps = 0.0;
  };

  /// @brief Centralized ground-segment ODTS *filter*, hosted on a `SurfaceStationManager`.
  ///
  /// A single EKF estimating every satellite's joint orbit + clock 8-state from
  /// `IslStationObs`s pushed by the `StationBeaconSensor`s on the physical surface-station
  /// agents -- no inter-satellite links. As onboard software it never synthesizes measurements
  /// from truth: it only consumes what the station sensors deliver and runs a predict/update
  /// cycle with its own (filter) dynamics. (Truth is read once, at setup, to seed the a-priori
  /// guess -- the analog of an uploaded a-priori ephemeris, as in the ground-station example.)
  /// Range-only ground tracking is weakly observable, so it uses a larger process noise and
  /// normalized-residual outlier rejection.
  class GroundOdtsApp : public Application {
  public:
    static constexpr int kSub = 8;

    GroundOdtsApp() = default;
    explicit GroundOdtsApp(Config& config);

    void Setup() override;
    void Step(Real t) override;
    void Log(Real /*t*/) override {}

    /// @brief Resolve a satellite name to its filter-state block index (for the sensors).
    int SatIndex(const std::string& sat_name) const;
    /// @brief Uniform-grid epoch index for time `t` [s].
    int EpochIndex(double t) const { return static_cast<int>(std::lround(t / dt_s_)); }
    /// @brief A station sensor pushes one observation into the current epoch's buffer.
    void AddMeasurement(const IslStationObs& m) { inbox_.push_back(m); }

    const std::vector<std::string>& SatelliteNames() const { return sat_names_; }
    const VecXd& TimeGrid() const { return t_grid_; }
    const MatXd& EstCentral() const { return est_central_; }  // [N x 8*n_sat]
    const std::vector<MatXd>& CovCentralFull() const {
      return cov_central_full_;
    }  // n_sat x [N x 64]

  protected:
    void Initialize();
    void RecordEpoch(int k);

    // Config
    int seed_ = 42;
    double dt_s_ = 60.0, duration_s_ = 21600.0;
    int moon_gravity_degree_filter_ = 8, moon_gravity_order_filter_ = 8;
    bool include_earth_ = true, include_sun_ = true, use_relativity_ = true;
    double integration_step_s_ = 60.0;
    double pseudorange_sigma_m_ = 1.0;
    bool include_station_doppler_ = true;
    double station_doppler_sigma_mps_ = 1.0e-3;
    double process_accel_sigma_mps2_ = 1.0e-6;
    double outlier_threshold_ = 3.0;
    double initial_position_sigma_m_ = 200.0, initial_velocity_sigma_mps_ = 0.1;
    double initial_clock_bias_sigma_s_ = 1.0e-6, initial_clock_drift_sigma_sps_ = 1.0e-9;
    std::vector<std::string> sat_names_, station_names_;

    // Runtime
    bool initialized_ = false;
    Real epoch0_ = 0.0;
    int n_sat_ = 0, n_state_ = 0;
    std::vector<Spacecraft*> sats_;         // for the a-priori seed only (not per-epoch)
    std::map<std::string, int> sat_index_;  // sat name -> filter block index
    std::vector<IslStationObs> inbox_;      // observations pushed by the station sensors
    Ptr<EKF> ekf_;
    std::mt19937 rng_;

    // Results
    VecXd t_grid_;
    MatXd est_central_;
    std::vector<MatXd> cov_central_full_;
  };

}  // namespace lupnt
