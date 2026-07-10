#pragma once

#include <string>
#include <vector>

#include "lupnt/applications/application.h"
#include "lupnt/core/config.h"
#include "lupnt/core/definitions.h"
#include "lupnt/numerics/filters/ekf.h"

namespace lupnt {

  class IslSatellite;

  /// @brief Centralized ground-segment ODTS filter, hosted on a `SurfaceStationManager`.
  ///
  /// A single EKF estimating every satellite's joint orbit + clock 8-state from the surface
  /// stations' one-way pseudorange (and optional Doppler) measurements alone -- no
  /// inter-satellite links. It pulls the satellites' and stations' truth via the agent graph,
  /// synthesizes the visible station pseudoranges each epoch, and runs a predict/update cycle.
  /// This is the ground-based baseline the distributed onboard filters (`SatelliteOdtsApp`) are
  /// compared against. Range-only ground tracking is weakly observable, so it uses a larger
  /// process noise and normalized-residual outlier rejection.
  class GroundOdtsApp : public Application {
  public:
    static constexpr int kSub = 8;

    GroundOdtsApp() = default;
    explicit GroundOdtsApp(Config& config);

    void Setup() override;
    void Step(Real t) override;
    void Log(Real /*t*/) override {}

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
    std::vector<IslSatellite*> sats_;
    std::vector<Vec3> stations_bf_;
    std::vector<double> station_mask_deg_;
    Ptr<EKF> ekf_;
    std::mt19937 rng_;

    // Results
    VecXd t_grid_;
    MatXd est_central_;
    std::vector<MatXd> cov_central_full_;
  };

}  // namespace lupnt
