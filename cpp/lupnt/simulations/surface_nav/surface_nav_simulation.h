#pragma once

#include <string>
#include <vector>

#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief Initial Cartesian state of one LCRNS relay satellite, defined at
  /// `SurfaceNavConfig::start_epoch_utc` in the Moon-centered inertial frame
  /// (`Frame::MOON_CI`), e.g. taken from the NASA LCRNS reference constellation.
  struct LcrnsSatConfig {
    std::string name = "SV";
    Vec3d r0_m = Vec3d::Zero();    ///< initial position [m], Frame::MOON_CI.
    Vec3d v0_mps = Vec3d::Zero();  ///< initial velocity [m/s], Frame::MOON_CI.
  };

  /// @brief Configuration for `RunSurfaceNav`: a lunar-surface rover fusing a full IMU
  /// (accelerometer + gyroscope, Kalibr noise model), LCRNS (LANS) pseudoranges, and a DEM
  /// altitude constraint, with the IMU biases estimated online.
  ///
  /// The rover truth path is defined in a local East-North-Up tangent plane centered on the
  /// selected DEM site (`site_lat_deg`/`site_lon_deg`); its Up coordinate follows the loaded
  /// DEM terrain. The body frame is x=forward (heading), y=left, z=up. LCRNS satellites are
  /// propagated as Keplerian orbits and their geometry is evaluated in the Moon-fixed frame.
  struct SurfaceNavConfig {
    int seed = 42;
    std::string start_epoch_utc = "2027-03-01T00:00:00";
    double duration_s = 1800.0;  ///< total arc [s].
    double dt_s = 1.0;           ///< IMU / filter / measurement step [s].

    // ---- Site / DEM -----------------------------------------------------------
    double site_lat_deg = -89.45;      ///< query latitude [deg] (default Site01).
    double site_lon_deg = 222.8;       ///< query east longitude [deg].
    double dem_half_width_m = 4000.0;  ///< half-width of the DEM crop window [m].
    double dem_max_res_m = 20.0;       ///< DEM downsample target spacing [m].

    // ---- Rover truth path (local ENU tangent plane) ---------------------------
    double rover_start_east_m = 0.0;         ///< start East offset from tile center [m].
    double rover_start_north_m = -250.0;     ///< start North offset from tile center [m].
    double rover_speed_mps = 2.0;            ///< horizontal ground speed [m/s].
    double rover_heading_deg = 0.0;          ///< initial heading (0=E, 90=N) [deg].
    double rover_turn_rate_dps = 0.6;        ///< heading rate [deg/s]; arcs sweep heading
                                             ///< (a full loop aids yaw / horizontal-bias obs).
    double rover_clock_bias_s = 1.0e-6;      ///< truth rover clock bias at t0 [s].
    double rover_clock_drift_sps = 1.0e-11;  ///< truth rover clock drift [s/s].

    // ---- IMU (Kalibr noise model) --------------------------------------------
    // https://github.com/ethz-asl/kalibr/wiki/IMU-Noise-Model. Defaults ~ LN200S-class.
    double accel_noise_density = 3.4e-4;  ///< accel white noise [m/s^2 / sqrt(Hz)].
    double accel_bias_rw = 1.0e-4;        ///< accel bias random walk [m/s^3 / sqrt(Hz)].
    double gyro_noise_density = 3.4e-6;   ///< gyro white noise [rad/s / sqrt(Hz)].
    double gyro_bias_rw = 1.0e-6;         ///< gyro bias random walk [rad/s^2 / sqrt(Hz)].
    double accel_bias0 = 1.0e-2;          ///< truth initial accel bias 1-sigma/axis [m/s^2].
    double gyro_bias0 = 5.0e-4;           ///< truth initial gyro bias 1-sigma/axis [rad/s].

    // ---- LCRNS constellation --------------------------------------------------
    std::vector<LcrnsSatConfig> satellites;  ///< relay satellites (e.g. 5).
    double elevation_mask_deg = 5.0;         ///< local-horizon mask for visibility [deg].
    double pseudorange_sigma_m = 1.0;        ///< receiver pseudorange noise 1-sigma [m].
    double sise_m = 3.0;                     ///< per-satellite signal-in-space error 1-sigma [m].

    // ---- Filter initialization & tuning --------------------------------------
    double init_pos_sigma_m = 50.0;
    double init_vel_sigma_mps = 0.1;
    double init_att_sigma_deg = 2.0;        ///< initial attitude uncertainty [deg].
    double init_accel_bias_sigma = 2.0e-2;  ///< initial accel-bias uncertainty [m/s^2].
    double init_gyro_bias_sigma = 1.0e-3;   ///< initial gyro-bias uncertainty [rad/s].
    double init_clock_bias_sigma_s = 1.0e-6;
    double init_clock_drift_sigma_sps = 1.0e-9;
    double dem_sigma_m = 5.0;  ///< DEM constraint pseudo-measurement 1-sigma [m].

    /// Scale applied to the filter's IMU bias random-walk densities relative to the truth
    /// values (>= 1). A modest inflation keeps the EKF consistent (error within 3-sigma) in
    /// the weakly-observable yaw / horizontal-accel-bias directions of a slow surface rover.
    double filter_bias_rw_scale = 3.0;

    bool enable_dem_constraint = true;
  };

  /// @brief Outputs of `RunSurfaceNav`: the loaded terrain plus per-epoch truth/estimate
  /// error and covariance (position, clock, attitude, and IMU biases), in a form directly
  /// plottable from Python (all matrices are plain `double`/`int`).
  struct SurfaceNavResults {
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

    // Time series (length N).
    VecXd time_s;
    MatXd pos_err_enu;       ///< N x 3, (truth - est) position error in ENU [m].
    MatXd pos_sigma_enu;     ///< N x 3, 1-sigma position uncertainty in ENU [m].
    VecXd pos_err_norm;      ///< N, 3D position error magnitude [m].
    VecXd clock_bias_err;    ///< N, clock-bias error [s].
    VecXd clock_bias_sigma;  ///< N, clock-bias 1-sigma [s].
    VecXi n_visible;         ///< N, number of visible LCRNS satellites.

    // IMU estimation (body frame), length N.
    MatXd accel_bias_err;    ///< N x 3, (truth - est) accel bias [m/s^2].
    MatXd accel_bias_sigma;  ///< N x 3, accel-bias 1-sigma [m/s^2].
    MatXd gyro_bias_err;     ///< N x 3, (truth - est) gyro bias [rad/s].
    MatXd gyro_bias_sigma;   ///< N x 3, gyro-bias 1-sigma [rad/s].
    MatXd att_err_deg;       ///< N x 3, attitude error (rotation vector) [deg].
    MatXd att_sigma_deg;     ///< N x 3, attitude 1-sigma [deg].

    MatXd rover_track_enu_truth;  ///< N x 2, truth (East, North) [m].
    MatXd rover_track_enu_est;    ///< N x 2, estimated (East, North) [m].
    VecXd rover_alt_truth;        ///< N, truth Up / elevation [m].

    std::vector<std::string> satellite_names;
  };

  /// @brief Run the surface-rover IMU + LCRNS + DEM navigation simulation end to end.
  ///
  /// Loads the DEM for the configured site (downloading from NASA PGDA if not cached),
  /// generates the rover truth trajectory and attitude over the terrain, propagates the
  /// LCRNS relays, synthesizes noisy accelerometer/gyroscope (Kalibr model) plus pseudorange
  /// measurements and the DEM altitude constraint, runs the `SurfaceRoverNavApp` strapdown INS
  /// EKF (estimating IMU biases online), and returns the logged truth/estimate/covariance
  /// series.
  SurfaceNavResults RunSurfaceNav(const SurfaceNavConfig& config);

}  // namespace lupnt
