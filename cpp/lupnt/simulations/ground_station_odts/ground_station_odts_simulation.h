#pragma once

#include <string>
#include <vector>

#include "lupnt/numerics/filters/batch_filter.h"
#include "lupnt/simulations/simulation.h"

namespace lupnt {

  /// @brief Geodetic location of a single ground-station antenna used by
  /// `GroundStationOdtsSimulation`. Converted to Earth-fixed (`Frame::ECEF`)
  /// Cartesian coordinates internally via `LatLonAltToCart`.
  struct GroundStationOdtsStationConfig {
    std::string name = "STATION";
    double latitude_deg = 0.0;
    double longitude_deg = 0.0;
    double altitude_m = 0.0;
  };

  /// @brief Default three-station Deep Space Network (DSN) 70 m antenna subnet --
  /// Goldstone (DSS-14), Canberra (DSS-43), Madrid (DSS-63) -- used as
  /// `GroundStationOdtsConfig::ground_stations` unless overridden.
  std::vector<GroundStationOdtsStationConfig> DefaultDsnStations();

  struct GroundStationOdtsConfig {
    int seed = 42;
    std::string start_epoch_utc = "2026-01-01T00:00:00";

    // Initial orbit: classical elements in the Moon-centered Orbital-Plane frame
    // (Frame::MOON_OP), converted to Frame::MOON_CI at start_epoch_utc.
    double orbit_a_m = 6541.4e3;
    double orbit_ecc = 0.6;
    double orbit_inc_rad = 56.2 * RAD;
    double orbit_raan_rad = 0.0;
    double orbit_argp_rad = 90.0 * RAD;
    double orbit_mean_anomaly_rad = 0.0;

    // Force model (the same dynamics are used for the truth trajectory and the
    // batch filter's design matrix; no truth/filter model mismatch).
    int moon_gravity_degree = 2;
    int moon_gravity_order = 2;
    bool include_earth = true;
    bool include_sun = true;
    bool use_relativity = false;
    bool use_srp = true;
    double srp_cr = 1.8;
    double srp_area_m2 = 0.02;
    double srp_mass_kg = 10.0;
    double integration_step_s = 10.0;

    // Ground stations and visibility.
    std::vector<GroundStationOdtsStationConfig> ground_stations = DefaultDsnStations();
    double elevation_mask_deg = 10.0;

    // Simulation timing.
    double duration_s = 36.0 * 3600.0;
    double obs_interval_s = 300.0;

    // Measurement types and noise (1-sigma).
    bool use_range = true;
    bool use_range_rate = true;
    double range_sigma_m = 10.0;
    double range_rate_sigma_mps = 1.0e-3;

    // A-priori state error injected into the batch filter's starting guess,
    // sampled ~ N(0, sigma^2) per component and seeded by `seed`.
    double initial_position_sigma_m = 2000.0;
    double initial_velocity_sigma_mps = 0.5;

    // Batch (iterative weighted-least-squares) filter -- see
    // lupnt/numerics/filters/batch_filter.h.
    int batch_max_iterations = 10;
    double batch_convergence_tol = 1.0e-3;
    bool batch_use_weights = true;
    bool batch_use_initialization = false;

    // Batch design-matrix Jacobian d(range, range_rate)/d(x0). When true (default),
    // the Jacobian is built analytically -- an autodiff state-transition matrix
    // Phi(t_i, t0) chained with the closed-form range/range-rate observation
    // partials -- so no finite differencing is used. When false, the legacy
    // finite-difference (numerical) Jacobian is used instead.
    bool batch_use_analytic_jacobian = true;

    // Square-root information filter (SRIF) + Dyer--McReynolds smoother, run after
    // the batch as a numerically-robust, process-noise-aware post-processor -- see
    // lupnt/numerics/filters/srif.h. Linearized about the converged batch
    // trajectory, it processes the same measurements once forward (sequential
    // square-root form) and smooths backward, yielding a per-epoch covariance that
    // reflects all measurements plus injected process noise -- a more realistic
    // (less optimistic) uncertainty than mapping the batch epoch covariance forward
    // with Phi P Phi^T.
    bool run_srif = true;
    // Continuous white-noise-acceleration process noise between epochs, modeling
    // unmodeled dynamics (higher-order gravity, SRP mismodeling, ...). With
    // srif_use_process_noise = false (or srif_accel_psd = 0) the SRIF is a pure
    // square-root batch and the smoother reproduces the propagated batch solution.
    bool srif_use_process_noise = true;
    double srif_accel_psd = 3.0e-13;  // per-axis PSD [m^2/s^3]
  };

  /// @brief Time series and batch-filter results of a `GroundStationOdtsSimulation`
  /// run, in SI units, `Frame::MOON_CI` unless noted otherwise.
  struct GroundStationOdtsResults {
    std::vector<std::string> station_names;

    // Truth trajectory over the full simulation grid.
    VecXd t_tdb;        // TDB epoch [s since J2000], size [N]
    MatXd truth_state;  // [N x 6], columns [r_x,r_y,r_z,v_x,v_y,v_z]

    // Visibility: topocentric elevation [deg] of the satellite as seen from each
    // station, evaluated on the same grid as t_tdb.
    MatXd elevation_deg;  // [N x n_stations]

    // Simulated tracking measurements (one row per station pass above the
    // elevation mask); obs_epoch_index indexes into t_tdb/truth_state/elevation_deg
    // and obs_station_index indexes into station_names.
    VecXi obs_epoch_index;
    VecXi obs_station_index;
    VecXd obs_range_true_m;
    VecXd obs_range_rate_true_mps;
    VecXd obs_range_m;         // range_true + noise (NaN if use_range == false)
    VecXd obs_range_rate_mps;  // range_rate_true + noise (NaN if use_range_rate == false)

    // Batch filter convergence history, one row per iteration performed.
    MatXd iteration_state_estimate;      // [K x 6], state at the *start* of each iteration
    VecXd iteration_correction_norm;     // [K]
    VecXd iteration_weighted_rms;        // [K], sqrt(mean(weight * residual^2)) over all obs
    VecXd iteration_pos_error_m;         // [K], iteration_state_estimate vs. x0_true
    VecXd iteration_vel_error_mps;       // [K]
    VecXd iteration_rms_range_m;         // [K]
    VecXd iteration_rms_range_rate_mps;  // [K]

    // Final batch-filter solution.
    Vec6d x0_true;
    Vec6d x0_initial_guess;
    Vec6d x0_estimated;
    MatXd covariance;  // [6 x 6] formal covariance of x0_estimated
    bool converged = false;
    int num_iterations = 0;

    // Estimated trajectory (x0_estimated propagated over the same grid as
    // truth_state), for direct comparison against truth_state.
    MatXd estimated_state;  // [N x 6]

    // Formal covariance of the estimated state propagated over the full grid,
    // P(t_i) = Phi(t_i, t0) * covariance * Phi(t_i, t0)^T, using the autodiff
    // state-transition matrix Phi. Each row is a row-major flattening of the
    // 6x6 covariance at that epoch (Frame::MOON_CI), so
    // estimated_covariance.row(i).reshaped(6, 6) recovers P(t_i).
    MatXd estimated_covariance;  // [N x 36]

    // Square-root information filter (SRIF) + smoother results over the full grid
    // (Frame::MOON_CI), populated when config.run_srif is true (empty otherwise).
    // The SRIF is linearized about the batch trajectory and processes the same
    // measurements once forward (filter) and backward (smoother) with injected
    // process noise. Covariance rows are row-major flattenings of the 6x6 block at
    // that epoch, matching estimated_covariance's layout.
    MatXd srif_filtered_state;       // [N x 6], forward-pass (filter) estimate
    MatXd srif_filtered_covariance;  // [N x 36], forward-pass covariance
    MatXd srif_smoothed_state;       // [N x 6], smoothed (all-data) estimate
    MatXd srif_smoothed_covariance;  // [N x 36], smoothed covariance
  };

  /// @brief Ground-station orbit determination and timing system (ODTS) simulation
  /// for a lunar satellite tracked by a network of Earth-based ground stations
  /// (e.g. the Deep Space Network) via two-way range and range-rate (Doppler).
  ///
  /// Propagates a truth trajectory (`NBodyDynamics`), runs a topocentric
  /// elevation-mask visibility analysis against each configured ground station,
  /// simulates noisy range/range-rate measurements over the visible passes, and
  /// recovers the trajectory from a deliberately perturbed initial guess with the
  /// iterative batch (weighted least-squares) filter in `lupnt/numerics/filters/batch_filter.h`.
  /// By default the batch design matrix is built analytically: the closed-form
  /// range/range-rate observation partials chained with the autodiff state-transition
  /// matrix `Phi(t_i, t0)` from `NBodyDynamics::Propagate(..., stm)` (no finite
  /// differencing). Setting `batch_use_analytic_jacobian = false` selects the legacy
  /// finite-difference Jacobian instead.
  class GroundStationOdtsSimulation : public Simulation {
  public:
    explicit GroundStationOdtsSimulation(GroundStationOdtsConfig config);

    void Setup() override;
    void Precompute() override;
    void Run() override;

    const GroundStationOdtsConfig& GetConfig() const { return config_; }
    const GroundStationOdtsResults& GetResults() const { return results_; }

  private:
    GroundStationOdtsConfig config_;
    GroundStationOdtsResults results_;
    bool setup_complete_ = false;
    bool precompute_complete_ = false;

    // Populated by Setup().
    Real t0_tdb_ = 0.0;
    VecX t_tdb_grid_;                   // [N], Real copy of results_.t_tdb for dynamics calls
    std::vector<Vec3> station_r_ecef_;  // [n_stations], Frame::ECEF, station position (fixed)
  };

}  // namespace lupnt
