#pragma once

#include <random>
#include <string>
#include <vector>

#include "lupnt/applications/application.h"
#include "lupnt/core/config.h"
#include "lupnt/core/definitions.h"
#include "lupnt/numerics/filters/ekf.h"

namespace lupnt {

  class Spacecraft;
  class NBodyDynamics;

  /// @brief Satellite-to-satellite angles-only orbit determination application, hosted on an
  /// OBSERVER `Spacecraft`.
  ///
  /// The observer measures the unit line-of-sight (bearing) direction to a TARGET
  /// `Spacecraft` whose orbit/ephemeris is treated as KNOWN, and runs an EKF that estimates
  /// its OWN orbit 6-state `[r(3), v(3)]` (no clock -- angles are clock-independent) from
  /// those bearings (`SatBearingMeasurement`). Because the target position is known, this is
  /// a well-posed *landmark-style* navigation problem that converges (the ill-posed
  /// angles-only target-range estimation is deliberately NOT attempted).
  ///
  /// The scenario is self-contained: both truth trajectories are propagated from the agents'
  /// truth dynamics (captured at setup) and the whole arc is filtered in a single end-of-arc
  /// solve, so a Monte-Carlo ensemble of independent measurement/initial-error draws can be
  /// run in one simulation. Per-epoch truth/estimate/1-sigma (first MC run) and the ensemble
  /// position-error statistics are exposed for plotting.
  class AnglesOdtsApp : public Application {
  public:
    static constexpr int kState = 6;  ///< observer orbit state `[r, v]`

    AnglesOdtsApp() = default;
    explicit AnglesOdtsApp(Config& config);

    void Setup() override;
    void Step(Real /*t*/) override {}  // work is done in the single end-of-arc Solve
    void Log(Real /*t*/) override {}

    // ---- Result accessors ----
    const VecXd& TimeGrid() const { return t_grid_; }
    const MatXd& TruthState() const { return truth_state_; }        // [N x 6] observer truth
    const MatXd& EstState() const { return est_state_; }            // [N x 6] estimate (MC run 0)
    const MatXd& SigmaState() const { return sigma_state_; }        // [N x 6] 1-sigma (MC run 0)
    const VecXd& PositionErrorRms() const { return pos_err_rms_; }  // [N] RMS over MC runs
    const VecXd& VelocityErrorRms() const { return vel_err_rms_; }  // [N] RMS over MC runs
    /// @brief [N x 19] convenience block: `[t, truth(6), est(6), sigma(6)]` (MC run 0).
    const MatXd& Trajectory() const { return trajectory_; }
    double FinalPositionErrorM() const { return final_pos_err_m_; }  // RMS-over-MC, last epoch
    double RmsPositionErrorM() const { return rms_pos_err_m_; }      // RMS over converged tail
    int MonteCarloRuns() const { return monte_carlo_runs_; }

  protected:
    void Solve();
    Ptr<NBodyDynamics> BuildFilterDynamics() const;

    // ---- Config (application block) ----
    std::string target_name_;
    int seed_ = 42;
    double dt_s_ = 60.0;
    double duration_s_ = 0.0;                              // 0 -> take the simulation duration
    double angle_sigma_rad_ = 5.0 * 4.84813681109536e-06;  // default 5 arcsec [rad]
    int moon_gravity_degree_filter_ = 8, moon_gravity_order_filter_ = 8;
    bool include_earth_ = true, include_sun_ = true, use_relativity_ = true;
    double integration_step_s_ = 60.0;
    double process_accel_sigma_mps2_ = 1.0e-9;
    double initial_position_sigma_m_ = 1000.0, initial_velocity_sigma_mps_ = 1.0;
    double outlier_threshold_ = 1.0e6;  // effectively off (benign geometry)
    int monte_carlo_runs_ = 1;

    // ---- Resolved / runtime ----
    Spacecraft* observer_ = nullptr;
    Spacecraft* target_ = nullptr;
    Real epoch0_ = 0.0;
    bool solved_ = false;

    // Results
    VecXd t_grid_;
    MatXd truth_state_, est_state_, sigma_state_, trajectory_;
    VecXd pos_err_rms_, vel_err_rms_;
    double final_pos_err_m_ = 0.0, rms_pos_err_m_ = 0.0;
  };

}  // namespace lupnt
