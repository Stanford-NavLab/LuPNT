#pragma once

#include <string>
#include <vector>

#include "lupnt/applications/application.h"
#include "lupnt/core/config.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  class NBodyDynamics;

  /// @brief One tracking observation reported by a `GroundStationTrackingApp` to
  /// the `GroundStationManagerApp`. Range and range-rate are frame-invariant, so
  /// the station reports its own state in the world (inertial) frame and the
  /// manager builds the design matrix without needing the station's geometry.
  struct StationMeasurement {
    double t = 0.0;                     // sim-relative epoch [s]
    int station_id = 0;                 // index into GroundStationManagerApp::StationNames()
    int epoch_index = 0;                // index into the manager's uniform epoch grid
    Vec6d station_mci = Vec6d::Zero();  // station state in the world frame at t (x0-independent)
    bool has_range = false;
    bool has_range_rate = false;
    double range = 0.0;             // noisy [m]
    double range_rate = 0.0;        // noisy [m/s]
    double range_sigma = 0.0;       // 1-sigma [m]     (measurement weight/noise)
    double range_rate_sigma = 0.0;  // 1-sigma [m/s]
  };

  /// @brief Centralized orbit-determination application hosted on a
  /// `GroundStationManager` agent.
  ///
  /// The manager is the estimator of the ground segment: `GroundStationTrackingApp`s
  /// running on the individual `GroundStation` agents push their visibility-gated
  /// range/range-rate observations here (`AddMeasurement`) as the arc plays out.
  /// At end of arc a single `Solve` runs the combined weighted-least-squares batch
  /// filter over *all* stations' measurements (analytic design matrix: closed-form
  /// range/range-rate partials chained with the target dynamics' autodiff STM), then
  /// a square-root information filter + Dyer--McReynolds smoother for a
  /// process-noise-aware, per-epoch covariance. The target dynamics come from the
  /// shared `World` (`World::MakeDynamics`), so truth and filter share one force
  /// model by construction.
  class GroundStationManagerApp : public Application {
  public:
    GroundStationManagerApp() = default;
    explicit GroundStationManagerApp(Config& config);

    /// @brief Resolve the target agent and shared dynamics, build the uniform
    /// epoch grid, seed the batch's perturbed initial guess, and schedule the
    /// single end-of-arc `Solve`.
    void Setup() override;
    void Step(Real /*t*/) override {}  // manager is measurement-driven, not stepped
    void Log(Real t) override;

    /// @brief Register a tracking station with the manager, returning its station
    /// id (index into `StationNames()`). Called from `GroundStationTrackingApp::Setup`.
    int RegisterStation(const std::string& station_name);

    /// @brief Aggregate one observation from a tracking station.
    void AddMeasurement(const StationMeasurement& m);

    /// @brief Map a sim-relative epoch [s] to its index in the uniform grid (nearest).
    int EpochIndex(double t) const;

    /// @brief Run the centralized batch + SRIF/smoother over all aggregated
    /// measurements. Scheduled once after the last tracking `Step`; idempotent.
    void Solve();

    // ---- Result accessors (valid after Solve) ----
    bool HasSolved() const { return solved_; }
    bool Converged() const { return converged_; }
    int NumIterations() const { return num_iterations_; }
    int NumMeasurements() const { return static_cast<int>(meas_.size()); }
    const std::vector<std::string>& StationNames() const { return station_names_; }

    const Vec6d& X0True() const { return x0_true_; }
    const Vec6d& X0InitialGuess() const { return x0_guess_; }
    const Vec6d& X0Estimated() const { return x0_est_; }
    const Mat6d& Covariance() const { return covariance_; }

    // Full-arc time series on the uniform epoch grid (world frame).
    const VecXd& TimeGrid() const { return t_grid_; }                 // [N] sim-relative [s]
    const MatXd& TruthState() const { return truth_state_; }          // [N x 6]
    const MatXd& EstimatedState() const { return estimated_state_; }  // [N x 6]
    const MatXd& EstimatedCovariance() const { return estimated_covariance_; }  // [N x 36]

    // Batch iteration history.
    const MatXd& IterationState() const { return iter_state_; }                  // [K x 6]
    const VecXd& IterationPosError() const { return iter_pos_err_; }             // [K]
    const VecXd& IterationVelError() const { return iter_vel_err_; }             // [K]
    const VecXd& IterationCorrectionNorm() const { return iter_corr_norm_; }     // [K]
    const VecXd& IterationWeightedRms() const { return iter_weighted_rms_; }     // [K]
    const VecXd& IterationRmsRange() const { return iter_rms_range_; }           // [K]
    const VecXd& IterationRmsRangeRate() const { return iter_rms_range_rate_; }  // [K]

    // SRIF forward filter + smoother (empty if run_srif is false).
    const MatXd& SrifFilteredState() const { return srif_filtered_state_; }     // [N x 6]
    const MatXd& SrifFilteredCovariance() const { return srif_filtered_cov_; }  // [N x 36]
    const MatXd& SrifSmoothedState() const { return srif_smoothed_state_; }     // [N x 6]
    const MatXd& SrifSmoothedCovariance() const { return srif_smoothed_cov_; }  // [N x 36]

  protected:
    // Config
    std::string target_name_;
    int seed_ = 42;
    double initial_position_sigma_m_ = 2000.0;
    double initial_velocity_sigma_mps_ = 0.5;
    int batch_max_iterations_ = 10;
    double batch_convergence_tol_ = 1.0e-3;
    double obs_interval_s_ = 300.0;  // uniform grid spacing for the time-series outputs
    bool run_srif_ = true;
    bool srif_use_process_noise_ = true;
    double srif_accel_psd_ = 3.0e-13;

    // Optional filter-dynamics overrides (else the shared world force model). `filter_dynamics`
    // applies to both estimators; `batch_dynamics` / `sequential_dynamics` override one each.
    Config filter_dyn_node_, batch_dyn_node_, srif_dyn_node_;
    bool has_filter_dyn_ = false, has_batch_dyn_ = false, has_srif_dyn_ = false;

    // Resolved at Setup
    Ptr<NBodyDynamics> dynamics_;        // world force model (default; also supplies GetFrame)
    Ptr<NBodyDynamics> batch_dynamics_;  // dynamics for the batch filter
    Ptr<NBodyDynamics> srif_dynamics_;   // dynamics for the sequential SRIF/smoother
    Real epoch0_ = 0.0;                  // absolute TDB epoch of sim time t = 0
    std::vector<std::string> station_names_;

    // Aggregated observations
    std::vector<StationMeasurement> meas_;

    // Uniform epoch grid [N] (sim-relative [s]) for all time-series outputs.
    VecXd t_grid_;

    // Results
    Vec6d x0_true_ = Vec6d::Zero();
    Vec6d x0_guess_ = Vec6d::Zero();
    Vec6d x0_est_ = Vec6d::Zero();
    Mat6d covariance_ = Mat6d::Zero();
    bool converged_ = false;
    bool solved_ = false;
    int num_iterations_ = 0;

    MatXd truth_state_, estimated_state_, estimated_covariance_;
    MatXd iter_state_;
    VecXd iter_pos_err_, iter_vel_err_, iter_corr_norm_, iter_weighted_rms_;
    VecXd iter_rms_range_, iter_rms_range_rate_;
    MatXd srif_filtered_state_, srif_filtered_cov_, srif_smoothed_state_, srif_smoothed_cov_;
  };

}  // namespace lupnt
