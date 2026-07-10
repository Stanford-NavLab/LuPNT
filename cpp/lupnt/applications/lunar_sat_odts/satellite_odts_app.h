#pragma once

#include <map>
#include <string>
#include <vector>

#include "lupnt/applications/application.h"
#include "lupnt/applications/ephemeris/lunanet_sat_app.h"
#include "lupnt/applications/lunar_sat_odts/isl_odts_app.h"
#include "lupnt/core/config.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  class IslSatellite;

  /// @brief One satellite's broadcast posterior (own 8-state mean + 8x8 covariance),
  /// published over the simulation's pub/sub bus each exchange epoch and consumed by the
  /// other satellites' `SatelliteOdtsApp`s.
  struct IslPosteriorMsg {
    std::string sender;
    double t = 0.0;  // sim-relative epoch [s] of this broadcast (for propagate-to-now on receipt)
    VecXd mean;      // [8]
    MatXd cov;       // [8 x 8]
  };

  /// @brief Onboard distributed ISL ODTS flight application, hosted on an `IslSatellite`.
  ///
  /// This is the *distributed* counterpart of the centralized `IslOdtsCoordinatorApp`: instead
  /// of one node running all N filters, each satellite runs its own. On each scheduled `Step`
  /// this app:
  ///   1. generates its OWN two-way crosslink measurements to every neighbour satellite by
  ///      pulling their truth (orbit + clock) via the agent graph
  ///      (`IslSatellite::GetTruthStateAt`),
  ///   2. generates its surface-station aiding pseudoranges/Doppler (pulling station truth),
  ///   3. runs its onboard Schmidt-EKF (`IslOdtsApp`), and
  ///   4. at `consider_exchange_interval_s`, broadcasts its own posterior over the pub/sub bus
  ///      (`Simulation::Publish`) and fuses the neighbours' most-recently-received broadcasts into
  ///      its consider blocks (Covariance-Intersection or naive overwrite). Because deliveries are
  ///      scheduled events, each satellite fuses the *previous* exchange's broadcasts -- a
  ///      deterministic one-interval latency, exactly like a real broadcast-ephemeris link.
  ///
  /// Measurement noise is drawn independently on each satellite (a two-way link is measured on
  /// both ends), so results differ from -- but are statistically consistent with -- the
  /// centralized coordinator.
  class SatelliteOdtsApp : public Application {
  public:
    static constexpr int kSub = 8;  // per-satellite state size [r,v,cb,cd]

    SatelliteOdtsApp() = default;
    explicit SatelliteOdtsApp(Config& config);

    void Setup() override;
    void Step(Real t) override;
    void Log(Real /*t*/) override {}

    // Result accessors (own filter, over the run; each [N x *]).
    const std::string& SatName() const { return sat_name_; }
    const std::vector<std::string>& NeighborNames() const { return neighbor_names_; }
    const VecXd& TimeGrid() const { return t_grid_; }
    const MatXd& TruthState() const { return truth_state_; }   // [N x 8] own truth
    const MatXd& OwnEstimate() const { return own_est_; }      // [N x 8] own estimate
    const MatXd& OwnCovDiag() const { return own_cov_diag_; }  // [N x 8] own cov diagonal
    const MatXd& OwnCovFull() const { return own_cov_full_; }  // [N x 64] own 8x8 cov row-major
    int NumAnchorsTotal() const { return total_anchors_; }     // diagnostic: total station PRs used

  protected:
    void Initialize();
    void RecordEpoch(int k);
    void OnPosterior(const IslPosteriorMsg& msg);  // pub/sub callback (buffers by sender)

    // --- Config (application block) ---
    int seed_ = 42;
    double dt_s_ = 60.0;
    double duration_s_ = 21600.0;
    int moon_gravity_degree_filter_ = 8, moon_gravity_order_filter_ = 8;
    bool include_earth_ = true, include_sun_ = true, use_relativity_ = true;
    double integration_step_s_ = 60.0;
    double range_sigma_m_ = 1.0, range_rate_sigma_mps_ = 1.0e-3;
    bool include_time_transfer_ = true, include_frequency_transfer_ = true;
    double time_transfer_sigma_m_ = 1.0, frequency_transfer_sigma_mps_ = 1.0e-3;
    double pseudorange_sigma_m_ = 1.0;
    bool include_station_doppler_ = true;
    double station_doppler_sigma_mps_ = 1.0e-3;
    double process_accel_sigma_mps2_ = 1.0e-8;
    double initial_position_sigma_m_ = 200.0, initial_velocity_sigma_mps_ = 0.1;
    double initial_clock_bias_sigma_s_ = 1.0e-6, initial_clock_drift_sigma_sps_ = 1.0e-9;
    double consider_position_sigma_m_ = 500.0, consider_velocity_sigma_mps_ = 0.05;
    double consider_clock_bias_sigma_s_ = 1.0e-6, consider_clock_drift_sigma_sps_ = 1.0e-9;
    double consider_exchange_interval_s_ = 600.0;
    bool exchange_use_covariance_intersection_ = false;
    double exchange_ci_weight_ = -1.0;
    std::vector<std::string> neighbor_names_;  // explicit; else auto-discover other IslSatellites
    std::vector<std::string> station_names_;

    // --- Resolved / runtime ---
    bool initialized_ = false;
    IslSatellite* self_ = nullptr;
    std::string sat_name_;
    std::vector<IslSatellite*> neighbors_;                   // block b (1..) -> neighbor
    std::vector<std::pair<Vec3, std::string>> stations_bf_;  // (MOON_PA pos, name); + mask below
    std::vector<double> station_mask_deg_;
    Ptr<IslOdtsApp> filter_;
    Ptr<LunaNetSatApp> host_;  // wraps `filter_` (its Setup builds the Schmidt-EKF)
    Ptr<JointOrbitClockDynamics> exchange_dyn_;  // propagate stale broadcasts to the current epoch
    Real epoch0_ = 0.0;
    int n_blk_ = 0;  // own + neighbors
    double next_exchange_s_ = 0.0;
    int total_anchors_ = 0;
    std::mt19937 rng_;
    std::map<std::string, IslPosteriorMsg> inbox_;  // latest received posterior per neighbour

    // Results
    VecXd t_grid_;
    MatXd truth_state_, own_est_, own_cov_diag_, own_cov_full_;
  };

}  // namespace lupnt
