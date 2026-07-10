#pragma once

#include <string>
#include <vector>

#include "lupnt/simulations/simulation.h"

namespace lupnt {

  /// @brief Initial Cartesian state and clock state for one satellite in an
  /// `IslOdtsSimulation`, defined at `IslOdtsConfig::start_epoch_utc` in the
  /// Moon-centered inertial frame (`Frame::MOON_CI`).
  struct IslOdtsSatelliteConfig {
    std::string name = "SAT";

    /// Initial Cartesian position [m] and velocity [m/s] in `Frame::MOON_CI`, e.g.
    /// taken directly from a published reference-constellation ephemeris (ICRF is
    /// treated as equivalent to `Frame::MOON_CI` for this purpose).
    Vec3d r0_m = Vec3d::Zero();
    Vec3d v0_mps = Vec3d::Zero();

    double clock_bias_s = 0.0;
    double clock_drift_sps = 0.0;
  };

  /// @brief Simplified inter-satellite-link (ISL) crosslink budget, used to report a
  /// representative carrier-to-noise density ratio alongside the range/Doppler
  /// measurements. Not used by the estimator; provided so link performance can be
  /// assessed without a full antenna-gain-pattern/transmitter model (see
  /// `lupnt::GnssMeasurement`/`lupnt::Antenna` for that level of fidelity).
  struct IslLinkBudgetConfig {
    bool enabled = true;
    double tx_power_dbw = 5.0;  // ~3 W crosslink transmit power
    double tx_gain_dbi = 20.0;  // representative crosslink antenna gain
    double rx_gain_dbi = 20.0;
    double frequency_hz = 25.5e9;  // representative Ka-band crosslink carrier
    double system_noise_temp_k = 500.0;
  };

  /// @brief Lunar-surface-station one-way pseudorange aiding for the constellation.
  ///
  /// When enabled, a single ground station fixed on the lunar surface (known position,
  /// via lunar rotation, and a known/disciplined clock treated as the timing reference)
  /// serves the satellites one at a time in a round-robin **rotation**: at each
  /// measurement epoch it points at the next satellite in the cycle that is above its
  /// elevation mask and provides that satellite a single one-way pseudorange.
  ///
  /// The pseudorange is modeled as the geometric station-to-satellite range plus the
  /// satellite's own clock bias (`rho = |r_sat - r_station| + c * b_sat + noise`).
  /// Because the station's inertial (`Frame::MOON_CI`) state is known -- it is fixed in
  /// the Moon-fixed frame and rotates with the Moon -- this one measurement gives the
  /// served satellite direct *absolute* position and clock-bias observability, unlike
  /// the two-way crosslinks which are clock-bias-free and only constrain relative
  /// geometry. Being a lunar-local link, its geometry is well conditioned (no ~4 deg
  /// Earth-GPS degeneracy), so no artificial noise inflation is required.
  struct IslSurfaceStationConfig {
    bool enabled = false;

    std::string name = "GS-Shackleton";

    // Fixed geodetic location on the Moon (Moon-fixed principal-axis frame). The
    // default is near the lunar south pole, the LCRNS reference-mission focus region.
    double latitude_deg = -89.9;
    double longitude_deg = 0.0;
    double altitude_m = 0.0;

    /// One-way pseudorange measurement noise 1-sigma [m].
    double pseudorange_sigma_m = 5.0;

    /// A satellite is only served when its topocentric elevation above the station's
    /// local horizon exceeds this mask [deg]; the rotation skips to the next visible
    /// satellite otherwise.
    double elevation_mask_deg = 5.0;
  };

  struct IslOdtsConfig {
    int seed = 42;
    std::string start_epoch_utc = "2027-03-01T00:00:00";
    double duration_s = 6.0 * 3600.0;
    double dt_s = 60.0;
    double integration_step_s = 60.0;

    /// Constellation of N >= 2 satellites, fully cross-linked. Every satellite runs
    /// its own onboard `SchmidtEKF` in parallel (a distributed navigation scheme):
    /// filter `j` estimates satellite `j`'s own 8-state and carries each of the other
    /// N-1 satellites as a Schmidt consider-state block, fed by two-way crosslinks to
    /// all of them. `satellites[0]` is only distinguished as the reference used for the
    /// crosslink-geometry reporting arrays.
    std::vector<IslOdtsSatelliteConfig> satellites;

    IslLinkBudgetConfig link_budget;

    /// Lunar surface station(s), each a beacon that exchanges a one-way pseudorange with
    /// every satellite above its elevation mask (satellite-side pseudoranges aid the
    /// onboard distributed filters; station-side pseudoranges feed the centralized ground
    /// filter, see `enable_centralized_ground_filter`). If `surface_stations` is left
    /// empty the single `surface_station` below is used when enabled; populate
    /// `surface_stations` to model a multi-station ground network.
    IslSurfaceStationConfig surface_station;
    std::vector<IslSurfaceStationConfig> surface_stations;

    /// Two-way inter-satellite time transfer: in addition to the (clock-free) two-way
    /// crosslink range, each link also measures the range-equivalent clock-bias difference
    /// between its endpoints, with noise `time_transfer_sigma_m` (comparable to the
    /// two-way ranging noise). This makes the constellation's *relative* clocks observable
    /// over the crosslink mesh; the surface station anchors the absolute time.
    bool enable_two_way_time_transfer = true;
    double time_transfer_sigma_m = 1.0;

    /// Two-way inter-satellite frequency transfer: the rate companion to the time
    /// transfer, measuring the range-rate-equivalent clock-*drift* difference between each
    /// link's endpoints (noise `frequency_transfer_sigma_mps`, comparable to the crosslink
    /// range-rate noise). Makes the constellation's relative clock drift observable.
    bool enable_two_way_frequency_transfer = true;
    double frequency_transfer_sigma_mps = 1.0e-3;

    /// One-way **Doppler** (pseudorange-rate) on the surface-station beacon links, in
    /// addition to the pseudorange: `u . (v_sat - v_station) + C * clock_drift_sat`. Adds
    /// velocity + clock-drift observability, which the range-only ground filter badly needs.
    bool enable_station_doppler = true;
    double station_doppler_sigma_mps = 1.0e-3;

    /// Also run a single **centralized** EKF (hosted on the ground) that estimates every
    /// satellite's full `[r, v, clock_bias, clock_drift]` state from the surface-station
    /// one-way pseudoranges *only* -- no inter-satellite-link measurements are downlinked.
    /// This is the classic ground-tracking baseline the distributed ISL scheme improves on.
    bool enable_centralized_ground_filter = true;

    /// Process-noise 1-sigma acceleration [m/s^2] for the centralized ground filter. Range-
    /// only ground tracking is weakly observable in cross-range, so this is typically much
    /// larger than the ISL filters' `process_accel_sigma_mps2` to keep the filter from
    /// becoming over-confident and diverging under the truth/filter force-model mismatch.
    double central_process_accel_sigma_mps2 = 1.0e-6;

    /// Normalized-residual outlier-rejection threshold [sigma] for the centralized ground
    /// filter. Range/Doppler tracking of the fast low-perilune passes produces occasional
    /// large linearized residuals that can diverge a naive EKF; rejecting residuals beyond
    /// this many sigma keeps it stable. Set very large (e.g. 1e12) to disable.
    double central_outlier_threshold = 3.0;

    /// Interval [s] at which the parallel filters exchange state: each satellite
    /// broadcasts its posterior own estimate (mean + covariance) and every other
    /// filter fuses it into the matching consider block, turning static neighbor
    /// priors into filter-improved priors. Set <= 0 to disable the exchange (each
    /// filter then keeps its independent, uncorrected consider states).
    double consider_exchange_interval_s = 600.0;

    /// How a broadcast neighbor estimate is fused into a filter's consider block.
    ///
    /// The parallel filters are *correlated* -- they share the same two-way crosslink
    /// observations, and every prior exchange has already mixed their estimates -- but
    /// that correlation is unknown and intractable to track on a fully-connected mesh
    /// (cycles cause "data incest"). Naively overwriting a consider block with a
    /// neighbor's posterior and zeroing its cross-covariance (the `false` setting)
    /// assumes independence, double-counts the shared measurements, and makes the
    /// reported covariance optimistic.
    ///
    /// When `true` (default), the exchange instead uses **Covariance Intersection**
    /// (Julier & Uhlmann): the fused information is a convex combination
    /// `Y_f = w * Y_own + (1-w) * Y_broadcast` with `w` chosen to minimize `tr(P_f)`.
    /// CI is guaranteed consistent for *any* (unknown) cross-correlation, so it removes
    /// both the exchange-induced and shared-measurement optimism at the price of some
    /// conservatism. NOTE: CI keeps *fusion* consistent only if each filter's own
    /// posterior is already consistent -- size `process_accel_sigma_mps2` to cover the
    /// truth/filter force-model mismatch as well.
    bool exchange_use_covariance_intersection = true;

    /// Fixed Covariance-Intersection weight `w` in `[0, 1]` for the neighbor fusion
    /// (weight on the *local* consider block; `1-w` on the broadcast). Set to a negative
    /// value (the default) to instead pick `w` per block by minimizing the trace of the
    /// fused block covariance -- the standard CI weight selection. Only used when
    /// `exchange_use_covariance_intersection` is true.
    double exchange_ci_weight = -1.0;

    // Force model
    int moon_gravity_degree_truth = 20;
    int moon_gravity_order_truth = 20;
    int moon_gravity_degree_filter = 8;
    int moon_gravity_order_filter = 8;
    bool include_earth = true;
    bool include_sun = true;
    bool use_relativity = true;

    // Measurement noise (two-way range / Doppler collapsed to the geometric
    // relative range and range-rate, see `isl_odts_simulation.cc`), applied to
    // every crosslink identically.
    double range_sigma_m = 1.0;
    double range_rate_sigma_mps = 1.0e-3;

    // Onboard Schmidt-EKF tuning. The hub filter estimates its own 8-state
    // [r, v, clock_bias, clock_drift] and considers each linked satellite's
    // 8-state broadcast/prior state without correcting it.
    double initial_position_sigma_m = 200.0;
    double initial_velocity_sigma_mps = 0.1;
    double initial_clock_bias_sigma_s = 1.0e-6;
    double initial_clock_drift_sigma_sps = 1.0e-9;

    // A satellite's knowledge of a neighbor's state -- e.g. from a periodically
    // refreshed broadcast/prior ephemeris -- is usually somewhat better than its
    // own just-initialized uncertainty but still not directly correctable onboard;
    // these sigmas seed the (fixed-mean) consider-state covariance block, applied
    // identically to every linked satellite.
    double consider_position_sigma_m = 100.0;
    double consider_velocity_sigma_mps = 0.05;
    double consider_clock_bias_sigma_s = 1.0e-7;
    double consider_clock_drift_sigma_sps = 1.0e-10;

    double process_accel_sigma_mps2 = 1.0e-8;
  };

  /// @brief Time series results of an `IslOdtsSimulation` run, in SI units,
  /// `Frame::MOON_CI`. Row `k` of every matrix corresponds to `t_s(k)`.
  struct IslOdtsResults {
    VecXd t_s;  // Elapsed time since start_epoch_utc [s], size [N]

    std::vector<std::string> satellite_names;  // size [n_sat]

    // Truth trajectories, one per satellite (same order as satellite_names), columns
    // [r_x,r_y,r_z,v_x,v_y,v_z,clock_bias_s,clock_drift_sps]
    std::vector<MatXd> truth_states;  // size [n_sat], each [N x 8]

    // Per-satellite onboard Schmidt-EKF estimate / covariance-diagonal. `est[j]` is
    // satellite j's own filter, stacked [own(8), consider blocks(8 each)] where own is
    // satellite j and the consider blocks are the other satellites in ascending global
    // index order (skipping j). Slice columns [0,8) for filter j's own-state estimate.
    std::vector<MatXd> est;       // size [n_sat], each [N x 8*n_sat]
    std::vector<MatXd> cov_diag;  // size [n_sat], each [N x 8*n_sat]

    // Full 8x8 own-state covariance of each satellite's own filter, row-major flattened
    // (64 columns per epoch, reshape to (N, 8, 8)). Unlike `cov_diag` this keeps the
    // off-diagonal terms, so it supports covariance-consistency (NEES) diagnostics that
    // need the whole own-block information matrix, not just its variances.
    std::vector<MatXd> cov_own_full;  // size [n_sat], each [N x 64]

    // Crosslink-geometry reporting (relative to satellites[0]): two-way range/Doppler
    // truth + noisy observation, one column per crosslink from satellite 0 to
    // satellite i+1 (n_links = n_sat - 1).
    MatXd range_true_m;         // [N x n_links]
    MatXd range_rate_true_mps;  // [N x n_links]
    MatXd range_obs_m;          // [N x n_links]
    MatXd range_rate_obs_mps;   // [N x n_links]

    // Pre-fit crosslink range residuals of each satellite's own filter, evaluated just
    // before that epoch's Update() -- the standard OD tracking diagnostic. `range_resid_m[j]`
    // holds filter j's residuals to its n_links neighbors (columns in ascending global
    // index order, skipping j).
    std::vector<MatXd> range_resid_m;  // size [n_sat], each [N x n_links]

    // Simplified crosslink carrier-to-noise density ratio [dB-Hz] (NaN if disabled),
    // relative to satellites[0] (same column convention as range_true_m).
    MatXd cn0_dbhz;  // [N x n_links]

    // Two-way inter-satellite time-transfer (range-equivalent clock-bias difference)
    // truth + noisy observation, relative to satellites[0] (same column convention as
    // range_true_m). All NaN if enable_two_way_time_transfer is false.
    MatXd time_transfer_true_m;  // [N x n_links]
    MatXd time_transfer_obs_m;   // [N x n_links]

    // --- Lunar surface station beacon(s). Each station exchanges a one-way pseudorange
    // with every satellite above its elevation mask (elevation-gated visibility).
    std::vector<MatXd> station_pos_mci;  // [n_station] each [N x 3] inertial pos [m], MOON_CI
    MatXd station_visible;               // [N x n_sat] # of stations seeing each satellite
    // Onboard distributed filter's station pseudorange (from the first station that sees
    // each satellite; NaN where no station sees it), per satellite:
    MatXd station_pr_true_m;   // [N x n_sat] geometric+clock truth pseudorange
    MatXd station_pr_obs_m;    // [N x n_sat] noisy pseudorange observation
    MatXd station_pr_resid_m;  // [N x n_sat] onboard filter's pre-fit residual

    // --- Centralized ground filter (populated only if enable_centralized_ground_filter):
    // a single EKF estimating all satellites from the station pseudoranges alone (no ISL).
    // est_central columns are [r,v,clock_bias,clock_drift] per satellite, in global order.
    MatXd est_central;                    // [N x 8*n_sat]
    std::vector<MatXd> cov_central_full;  // [n_sat] each [N x 64] own 8x8 cov (row-major)
  };

  // The distributed ISL ODTS run is driven by `IslOdtsCoordinatorApp`
  // (applications/lunar_sat_odts/isl_odts_coordinator_app.h) on a thin `IslOdtsManager`
  // agent, via `pnt.Simulation(...)`. The `IslOdtsConfig`/`IslOdtsResults` structs above
  // are the shared config/result payloads consumed by that app.

}  // namespace lupnt
