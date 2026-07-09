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
  class Dynamics;

  /// @brief Ground-station orbit-determination application.
  ///
  /// Runs on a `GroundStation` agent and estimates the orbit of a *target* satellite
  /// agent (referenced by name) from two-way range and range-rate (Doppler) tracking.
  /// On each scheduled `Step` the app queries the target agent's truth state, checks
  /// topocentric visibility against an elevation mask, and — if visible — simulates a
  /// noisy range/range-rate measurement and stores it. After the tracking arc it runs
  /// the shared iterative batch (weighted least-squares) filter
  /// (`lupnt::RunBatchFilter`) with an *analytic* design matrix: the closed-form
  /// range/range-rate observation partials chained with the target dynamics' autodiff
  /// state-transition matrix `Phi(t_i, t0)` (no finite differencing).
  ///
  /// This is the agent-based counterpart of `GroundStationOdtsSimulation`: the truth
  /// trajectory and station geometry come from first-class agents rather than being
  /// propagated internally, so multiple ground-station apps and/or target satellites
  /// can participate in a single `Simulation`.
  class GroundStationOdtsApp : public Application {
  public:
    GroundStationOdtsApp() = default;

    /// @brief Construct from the `application:` YAML block of a `GroundStation` agent.
    /// Recognized keys (besides the base `name`/`frequency`): `target` (target agent
    /// name), `elevation_mask_deg`, `use_range`, `use_range_rate`, `range_sigma_m`,
    /// `range_rate_sigma_mps`, `seed`, `initial_position_sigma_m`,
    /// `initial_velocity_sigma_mps`, `batch_max_iterations`, `batch_convergence_tol`.
    explicit GroundStationOdtsApp(Config& config);

    /// @brief Resolve the target agent and station geometry, seed the batch filter's
    /// perturbed initial guess, schedule periodic measurement `Step`s (via the base
    /// class) and a single end-of-arc batch `Solve`.
    void Setup() override;

    /// @brief One measurement epoch: query the target agent, gate on the elevation
    /// mask, and (if visible) store a noisy range/range-rate measurement.
    void Step(Real t) override;

    /// @brief Log the batch-filter estimate/covariance once solved.
    void Log(Real t) override;

    /// @brief Run the batch orbit-determination over all accumulated measurements.
    /// Scheduled to run once after the last measurement `Step`; idempotent.
    void Solve();

    // Result accessors (valid after Solve()).
    bool HasSolved() const { return solved_; }
    bool Converged() const { return converged_; }
    int NumIterations() const { return num_iterations_; }
    int NumMeasurements() const { return static_cast<int>(meas_t_.size()); }
    const Vec6d& X0True() const { return x0_true_; }
    const Vec6d& X0InitialGuess() const { return x0_guess_; }
    const Vec6d& X0Estimated() const { return x0_est_; }
    const MatXd& Covariance() const { return covariance_; }

  protected:
    // Config
    std::string target_name_;
    double elevation_mask_deg_ = 10.0;
    bool use_range_ = true;
    bool use_range_rate_ = true;
    double range_sigma_m_ = 10.0;
    double range_rate_sigma_mps_ = 1.0e-3;
    int seed_ = 42;
    double initial_position_sigma_m_ = 2000.0;
    double initial_velocity_sigma_mps_ = 0.5;
    int batch_max_iterations_ = 10;
    double batch_convergence_tol_ = 1.0e-3;

    // Resolved at Setup()
    AgentWithDynamics* target_ = nullptr;  // tracked satellite agent
    Dynamics* target_dynamics_ = nullptr;  // its dynamics (autodiff STM source)
    Real epoch0_ = 0.0;                    // absolute TDB epoch of sim time t = 0
    Vec3 station_r_ = Vec3::Zero();        // station position, station body-fixed frame
    Frame station_frame_ = Frame::UNDEFINED;

    // Target truth state propagated incrementally with the target's dynamics as Step()
    // advances (so we never re-propagate the whole arc from t = 0 per measurement).
    Real target_t_ = 0.0;
    Vec6 target_rv_ = Vec6::Zero();  // MOON_CI truth state at target_t_

    // Accumulated measurements (sim-relative epoch [s], noisy values in SI)
    std::vector<double> meas_t_;
    std::vector<double> meas_range_;
    std::vector<double> meas_range_rate_;
    std::mt19937 noise_rng_;

    // Results
    Vec6d x0_true_ = Vec6d::Zero();
    Vec6d x0_guess_ = Vec6d::Zero();
    Vec6d x0_est_ = Vec6d::Zero();
    MatXd covariance_;
    bool converged_ = false;
    bool solved_ = false;
    int num_iterations_ = 0;
  };

}  // namespace lupnt
