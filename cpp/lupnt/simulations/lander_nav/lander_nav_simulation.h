#pragma once

#include <string>
#include <vector>

#include "lupnt/core/definitions.h"
#include "lupnt/simulations/surface_nav/surface_nav_simulation.h"  // LcrnsSatConfig

namespace lupnt {

  /// @brief Configuration for `RunLanderNav`: a lunar lander on powered descent to a south-pole
  /// site, navigating with an error-state INS EKF that fuses a full IMU (accelerometer +
  /// gyroscope, Kalibr noise model), a nadir radar altimeter (height above DEM terrain),
  /// crater-landmark bearings (terrain-relative navigation against a synthetic crater map), and
  /// LunaNet (LANS) pseudoranges, with the IMU biases estimated online.
  ///
  /// The descent truth path is defined in a local East-North-Up tangent plane centred on the
  /// selected DEM site: the horizontal position glides (via a smoothstep) from
  /// `(start_east, start_north)` to `(end_east, end_north)` while the height above terrain
  /// descends from `start_alt_m` to `end_alt_m`. The body frame is x=forward (heading), y=left,
  /// z=up; the down-looking camera boresight is nadir (body -z). LunaNet relays are propagated
  /// as Keplerian orbits and evaluated in the Moon-fixed frame.
  struct LanderNavConfig {
    int seed = 42;
    std::string start_epoch_utc = "2027-03-01T00:00:00";
    double duration_s = 300.0;  ///< total descent arc [s].
    double dt_s = 0.5;          ///< IMU / filter / measurement step [s].

    // ---- Site / DEM -----------------------------------------------------------
    double site_lat_deg = -89.45;      ///< query latitude [deg] (default Site01).
    double site_lon_deg = 222.8;       ///< query east longitude [deg].
    double dem_half_width_m = 4000.0;  ///< half-width of the DEM crop window [m].
    double dem_max_res_m = 20.0;       ///< DEM downsample target spacing [m].

    // ---- Lander descent truth (local ENU tangent plane) -----------------------
    double descent_start_east_m = -2500.0;    ///< start East offset (downrange) [m].
    double descent_start_north_m = 400.0;     ///< start North offset (crossrange) [m].
    double descent_end_east_m = 0.0;          ///< landing East offset [m].
    double descent_end_north_m = 0.0;         ///< landing North offset [m].
    double descent_start_alt_m = 2000.0;      ///< height above terrain at t0 [m].
    double descent_end_alt_m = 15.0;          ///< height above terrain at touchdown/hover [m].
    double descent_heading_deg = 0.0;         ///< body-x heading in ENU (0=E, 90=N) [deg].
    double lander_clock_bias_s = 1.0e-6;      ///< truth lander clock bias at t0 [s].
    double lander_clock_drift_sps = 1.0e-11;  ///< truth lander clock drift [s/s].

    /// Optional externally-supplied reference (truth) trajectory, as `N x 3` rows of local
    /// East-North-Up position [m] about the DEM tile center (U is height above the site datum,
    /// NOT above the terrain). When non-empty it **overrides** the built-in smoothstep descent
    /// (the `descent_*` fields above are ignored) and defines the epoch count `N` (so its rows
    /// must equal `round(duration_s / dt_s) + 1`). Populate it from any guidance law — e.g. the
    /// `pylupnt.lander_guidance` ZEM/ZEV or convex-optimization generators. Velocity, specific
    /// force, and angular rate are still derived by finite-differencing this path.
    MatXd ref_traj_enu;

    // ---- IMU (Kalibr noise model) --------------------------------------------
    // https://github.com/ethz-asl/kalibr/wiki/IMU-Noise-Model. Defaults ~ LN200S-class.
    double accel_noise_density = 3.4e-4;  ///< accel white noise [m/s^2 / sqrt(Hz)].
    double accel_bias_rw = 1.0e-4;        ///< accel bias random walk [m/s^3 / sqrt(Hz)].
    double gyro_noise_density = 3.4e-6;   ///< gyro white noise [rad/s / sqrt(Hz)].
    double gyro_bias_rw = 1.0e-6;         ///< gyro bias random walk [rad/s^2 / sqrt(Hz)].
    double accel_bias0 = 1.0e-2;          ///< truth initial accel bias 1-sigma/axis [m/s^2].
    double gyro_bias0 = 5.0e-4;           ///< truth initial gyro bias 1-sigma/axis [rad/s].

    // ---- Radar altimeter ------------------------------------------------------
    bool enable_altimeter = true;
    double altimeter_sigma_m = 2.0;         ///< altimeter noise 1-sigma [m].
    double altimeter_max_range_m = 3000.0;  ///< max altitude with a valid return [m].

    // ---- Crater-landmark camera (terrain-relative navigation) -----------------
    bool enable_craters = true;
    int n_craters = 60;                     ///< number of synthetic craters in the map.
    double crater_field_radius_m = 3000.0;  ///< craters scattered within this radius of the site.
    double camera_fov_deg = 30.0;           ///< camera half-cone about nadir [deg].
    int max_craters_per_epoch = 6;          ///< cap on landmarks tracked per epoch.
    double crater_sigma_arcsec = 120.0;     ///< bearing noise 1-sigma [arcsec] (~0.58 mrad).

    // ---- LunaNet (LCRNS) constellation ----------------------------------------
    bool enable_lunanet = true;
    std::vector<LcrnsSatConfig> satellites;  ///< relay satellites (e.g. 5).
    double elevation_mask_deg = 5.0;         ///< local-horizon mask for visibility [deg].
    double pseudorange_sigma_m = 1.0;        ///< receiver pseudorange noise 1-sigma [m].
    double sise_m = 3.0;                     ///< per-satellite signal-in-space error 1-sigma [m].

    // ---- Filter initialization & tuning --------------------------------------
    double init_pos_sigma_m = 100.0;
    double init_vel_sigma_mps = 1.0;
    double init_att_sigma_deg = 2.0;        ///< initial attitude uncertainty [deg].
    double init_accel_bias_sigma = 2.0e-2;  ///< initial accel-bias uncertainty [m/s^2].
    double init_gyro_bias_sigma = 1.0e-3;   ///< initial gyro-bias uncertainty [rad/s].
    double init_clock_bias_sigma_s = 1.0e-6;
    double init_clock_drift_sigma_sps = 1.0e-9;

    /// Scale applied to the filter's IMU bias random-walk densities relative to the truth
    /// values (>= 1), for EKF consistency in weakly-observable directions.
    double filter_bias_rw_scale = 3.0;
  };

  /// @brief Outputs of `RunLanderNav`: the loaded terrain and crater map plus per-epoch
  /// truth/estimate error and covariance (position, velocity, clock, attitude, IMU biases),
  /// in a form directly plottable from Python (all matrices are plain `double`/`int`).
  struct LanderNavResults {
    std::string site_id;
    std::string site_name;
    double site_lat_deg = 0.0;
    double site_lon_deg = 0.0;

    // Terrain grid (native DEM projected meters), for plotting.
    MatXd dem_x;
    MatXd dem_y;
    MatXd dem_elevation;
    double dem_center_x = 0.0;  ///< native x of the tile center (== ENU East origin).
    double dem_center_y = 0.0;  ///< native y of the tile center (== ENU North origin).

    // Synthetic crater map (ENU about tile center), for plotting.
    MatXd crater_enu;  ///< M x 2, crater (East, North) [m].

    // Time series (length N).
    VecXd time_s;
    MatXd pos_err_enu;       ///< N x 3, (truth - est) position error in ENU [m].
    MatXd pos_sigma_enu;     ///< N x 3, 1-sigma position uncertainty in ENU [m].
    VecXd pos_err_norm;      ///< N, 3D position error magnitude [m].
    VecXd vel_err_norm;      ///< N, 3D velocity error magnitude [m/s].
    VecXd clock_bias_err;    ///< N, clock-bias error [s].
    VecXd clock_bias_sigma;  ///< N, clock-bias 1-sigma [s].
    VecXi n_visible_sat;     ///< N, number of visible LunaNet satellites.
    VecXi n_craters;         ///< N, number of tracked crater landmarks.

    // IMU estimation (body frame), length N.
    MatXd accel_bias_err;    ///< N x 3, (truth - est) accel bias [m/s^2].
    MatXd accel_bias_sigma;  ///< N x 3, accel-bias 1-sigma [m/s^2].
    MatXd gyro_bias_err;     ///< N x 3, (truth - est) gyro bias [rad/s].
    MatXd gyro_bias_sigma;   ///< N x 3, gyro-bias 1-sigma [rad/s].
    MatXd att_err_deg;       ///< N x 3, attitude error (rotation vector) [deg].
    MatXd att_sigma_deg;     ///< N x 3, attitude 1-sigma [deg].
    MatXd att_quat_est;      ///< N x 4, estimated body-to-nav attitude quaternion [w,x,y,z] (MEKF).

    // Descent trajectory (ENU about tile center), length N.
    MatXd traj_enu_truth;  ///< N x 3, truth (East, North, Up) [m].
    MatXd traj_enu_est;    ///< N x 3, estimated (East, North, Up) [m].
    VecXd alt_truth;       ///< N, truth height above terrain [m].
    VecXd alt_est;         ///< N, estimated height above terrain [m].

    std::vector<std::string> satellite_names;
  };

  /// @brief Run the lunar-lander IMU + altimeter + crater-bearing + LunaNet navigation
  /// simulation end to end.
  ///
  /// Loads the DEM for the configured site (downloading from NASA PGDA if not cached), builds a
  /// synthetic crater-landmark map on the terrain, generates the powered-descent truth
  /// trajectory and attitude, propagates the LunaNet relays, then constructs a `Lander` agent
  /// hosting a `LanderNavApp` (attached via `Agent::SetApplication`) and drives its error-state
  /// INS EKF epoch by epoch, returning the logged truth/estimate/covariance series.
  LanderNavResults RunLanderNav(const LanderNavConfig& config);

}  // namespace lupnt
