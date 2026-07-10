#pragma once

#include <random>
#include <string>
#include <vector>

#include "lupnt/applications/application.h"
#include "lupnt/core/config.h"
#include "lupnt/core/definitions.h"
#include "lupnt/measurements/surface_measurements.h"

namespace lupnt {

  /// @brief Onboard tuning for `SurfaceRoverNavApp` (Kalibr IMU noise densities + measurement/
  /// clock process noise).
  ///
  /// The IMU noise follows the Kalibr model
  /// (https://github.com/ethz-asl/kalibr/wiki/IMU-Noise-Model): each sensor has a
  /// continuous-time white-noise density and a bias random-walk density. Defaults are
  /// representative tactical-grade (LN200S-like) values.
  struct SurfaceRoverNavAppParams {
    // ---- IMU (Kalibr) noise densities ----------------------------------------
    double accel_noise_density = 3.4e-4;  ///< accel white noise [m/s^2 / sqrt(Hz)]
    double accel_bias_rw = 1.0e-4;        ///< accel bias random walk [m/s^3 / sqrt(Hz)]
    double gyro_noise_density = 3.4e-6;   ///< gyro white noise [rad/s / sqrt(Hz)]
    double gyro_bias_rw = 1.0e-6;         ///< gyro bias random walk [rad/s^2 / sqrt(Hz)]

    // ---- Aiding measurement noise --------------------------------------------
    double pseudorange_sigma_m = 1.0;  ///< LANS receiver pseudorange noise 1-sigma [m]
    double sise_m = 3.0;               ///< LANS signal-in-space error 1-sigma [m]
    double dem_sigma_m = 5.0;          ///< DEM altitude-constraint 1-sigma [m]

    // ---- Clock process noise --------------------------------------------------
    double clock_bias_process_sigma = 1.0e-11;   ///< clock-bias random walk [s / sqrt(s)]
    double clock_drift_process_sigma = 1.0e-13;  ///< clock-drift random walk [s/s / sqrt(s)]
  };

  /// @brief Scenario configuration for the surface-rover navigation application: a lunar-surface
  /// rover fusing a full IMU (accelerometer + gyroscope, Kalibr noise model), LCRNS (LANS)
  /// pseudoranges, and a DEM altitude constraint, with the IMU biases estimated online.
  ///
  /// The rover truth path is defined in the local East-North-Up tangent plane centered on the
  /// `World`'s DEM site; its Up coordinate follows the loaded DEM terrain. The body frame is
  /// x=forward (heading), y=left, z=up. LCRNS satellites are propagated as Keplerian orbits
  /// app-internally and their geometry is evaluated in the Moon-fixed frame. The site / DEM
  /// itself is owned by the shared `World` (`world:` config block), not by this struct.
  struct SurfaceNavConfig {
    int seed = 42;
    std::string start_epoch_utc = "2027-03-01T00:00:00";
    double duration_s = 1800.0;  ///< total arc [s].
    double dt_s = 1.0;           ///< IMU / filter / measurement step [s].

    // ---- Site / DEM (site/half_width used only by the legacy free function; the agent-based
    // app takes the DEM from the shared `World`. `dem_max_res_m` is also the terrain-slope
    // finite-difference step and must match the World's `dem.max_res_m`). ------------------
    double site_lat_deg = -89.45;      ///< query latitude [deg] (default Site01).
    double site_lon_deg = 222.8;       ///< query east longitude [deg].
    double dem_half_width_m = 4000.0;  ///< half-width of the DEM crop window [m].
    double dem_max_res_m = 20.0;       ///< DEM downsample target spacing / slope FD step [m].

    // ---- Rover truth path (local ENU tangent plane) ---------------------------
    double rover_start_east_m = 0.0;         ///< start East offset from tile center [m].
    double rover_start_north_m = -250.0;     ///< start North offset from tile center [m].
    double rover_speed_mps = 2.0;            ///< horizontal ground speed [m/s].
    double rover_heading_deg = 0.0;          ///< initial heading (0=E, 90=N) [deg].
    double rover_turn_rate_dps = 0.6;        ///< heading rate [deg/s].
    double rover_clock_bias_s = 1.0e-6;      ///< truth rover clock bias at t0 [s].
    double rover_clock_drift_sps = 1.0e-11;  ///< truth rover clock drift [s/s].

    // ---- IMU (Kalibr noise model) --------------------------------------------
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
    /// values (>= 1). A modest inflation keeps the EKF consistent in weakly-observable
    /// yaw / horizontal-accel-bias directions of a slow surface rover.
    double filter_bias_rw_scale = 3.0;

    bool enable_dem_constraint = true;
  };

  /// @brief Per-epoch truth/estimate error and covariance series produced by the rover
  /// navigation application, in a form directly plottable from Python (all matrices are plain
  /// `double`/`int`). Mirrors the fields the old `RunSurfaceNav` returned.
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

  /// @brief Strapdown inertial-navigation error-state EKF that fuses a full IMU
  /// (accelerometer + gyroscope), LCRNS (LANS) pseudoranges, and a DEM altitude constraint
  /// to position and orient a surface rover, **estimating the IMU biases online**.
  ///
  /// The nominal state is `[r(3), v(3), R_b2n(SO3), b_a(3), b_g(3), clock_bias, clock_drift]`
  /// in a Moon-fixed frame; the 17-element error state
  /// `[dr, dv, dtheta, db_a, db_g, d(cb), d(cd)]` (see `kSurfaceNavErrorStateSize`) carries
  /// the covariance. Each step it `Predict`s with an IMU sample (bias-corrected strapdown
  /// mechanization + Kalibr process noise), then applies scalar updates: one per visible LANS
  /// satellite (`UpdateLans`) and, optionally, a linearized DEM altitude constraint
  /// (`UpdateScalar`).
  ///
  /// **Hosted on a `Rover` agent as its `Application`** (attach via `Agent::SetApplication`).
  /// When constructed from a YAML `application:` block, the app is self-driving: `Setup()`
  /// (called by the `Simulation` once the agent's `World` is available) precomputes the truth
  /// trajectory, IMU truth, relay orbits, and the perturbed initial estimate, seeds the filter,
  /// and schedules `Step(t)` at the IMU cadence; each `Step` synthesizes that epoch's IMU and
  /// pseudorange measurements from the shared `World` (terrain, gravity, ENU frame) and runs the
  /// predict/update cycle, recording the truth/estimate/covariance series accessible from the
  /// result accessors below. The lower-level `Configure`/`Predict`/`Update*` core can also be
  /// driven directly (e.g. from a custom driver) via the params-struct constructor.
  class SurfaceRoverNavApp : public Application {
  public:
    SurfaceRoverNavApp() = default;
    explicit SurfaceRoverNavApp(const SurfaceRoverNavAppParams& params) : params_(params) {}

    /// @brief Construct a self-driving app from the `application:` block of a `Rover` agent.
    /// Reads the scenario (`SurfaceNavConfig`) and IMU/measurement tuning; the site / DEM is
    /// taken from the shared `World` in `Setup()`.
    explicit SurfaceRoverNavApp(Config& config);

    /// @brief Precompute the truth trajectory, IMU truth, relay orbits, and perturbed initial
    /// estimate (in the same RNG draw order as the legacy `RunSurfaceNav`), seed the filter, log
    /// epoch 0, and schedule periodic `Step`s at the IMU cadence. Only active when the app was
    /// constructed from a YAML config; the source of terrain / gravity / ENU is `GetWorld()`.
    void Setup() override;

    /// @brief One self-driving epoch: set the host agent's truth state, synthesize + `Predict`
    /// the IMU sample, apply LANS pseudorange and DEM-constraint updates (same order as the
    /// legacy loop body), and record the result series. No-op for the params-struct driver.
    void Step(Real t) override;

    void Log(Real t) override;

    /// @brief Seed the filter nominal state and covariance.
    void Configure(double t0, const Vec3d& r0, const Vec3d& v0, const Mat3d& R0, const Vec3d& ba0,
                   const Vec3d& bg0, double cb0, double cd0, const MatXd& P0);

    /// @brief Predict to `t = t_ + dt` using a full IMU sample.
    void Predict(const SurfaceImuMeasurement& imu, double dt);

    /// @brief Apply a sequential scalar EKF update for each LANS pseudorange.
    void UpdateLans(const std::vector<SurfaceLansMeasurement>& meas);

    /// @brief Generic scalar EKF update in the error state, immediately injecting the resulting
    /// correction into the nominal state (error-state reset).
    void UpdateScalar(const VecXd& H, double z_pred, double z_meas, double variance);

    double time() const { return t_; }                ///< current filter epoch [s].
    const Vec3d& position() const { return r_; }      ///< estimated position [m].
    const Vec3d& velocity() const { return v_; }      ///< estimated velocity [m/s].
    const Mat3d& attitude() const { return R_b2n_; }  ///< estimated body-to-nav rotation.
    const Vec3d& accel_bias() const { return ba_; }   ///< estimated accel bias [m/s^2].
    const Vec3d& gyro_bias() const { return bg_; }    ///< estimated gyro bias [rad/s].
    double clock_bias() const { return cb_; }         ///< estimated clock bias [s].
    double clock_drift() const { return cd_; }        ///< estimated clock drift [s/s].
    const MatXd& covariance() const { return P_; }    ///< error-state covariance.
    const SurfaceRoverNavAppParams& params() const { return params_; }

    // ---- Result series (valid after the Simulation has run) -------------------
    const SurfaceNavConfig& config() const { return cfg_; }
    const SurfaceNavResults& results() const { return res_; }
    const std::string& site_id() const { return res_.site_id; }
    const std::string& site_name() const { return res_.site_name; }
    const MatXd& dem_x() const { return res_.dem_x; }
    const MatXd& dem_y() const { return res_.dem_y; }
    const MatXd& dem_elevation() const { return res_.dem_elevation; }
    const VecXd& time_series() const { return res_.time_s; }
    const MatXd& pos_err_enu() const { return res_.pos_err_enu; }
    const MatXd& pos_sigma_enu() const { return res_.pos_sigma_enu; }
    const VecXd& pos_err_norm() const { return res_.pos_err_norm; }
    const VecXd& clock_bias_err() const { return res_.clock_bias_err; }
    const VecXd& clock_bias_sigma() const { return res_.clock_bias_sigma; }
    const VecXi& n_visible() const { return res_.n_visible; }
    const MatXd& accel_bias_err() const { return res_.accel_bias_err; }
    const MatXd& accel_bias_sigma() const { return res_.accel_bias_sigma; }
    const MatXd& gyro_bias_err() const { return res_.gyro_bias_err; }
    const MatXd& gyro_bias_sigma() const { return res_.gyro_bias_sigma; }
    const MatXd& att_err_deg() const { return res_.att_err_deg; }
    const MatXd& att_sigma_deg() const { return res_.att_sigma_deg; }
    const MatXd& rover_track_enu_truth() const { return res_.rover_track_enu_truth; }
    const MatXd& rover_track_enu_est() const { return res_.rover_track_enu_est; }
    const VecXd& rover_alt_truth() const { return res_.rover_alt_truth; }
    const std::vector<std::string>& satellite_names() const { return res_.satellite_names; }

  private:
    /// Build the discrete process-noise covariance for a step of length `dt`.
    MatXd ProcessNoise(double dt) const;
    /// Inject a full error-state correction into the nominal state.
    void InjectErrorState(const VecXd& dx);
    /// Precompute the truth trajectory, IMU truth, relay orbits, and perturbed initial estimate;
    /// seed the filter and log epoch 0. Deferred to the first `Step` so any programmatic config
    /// set between `Simulation` construction and `run()` is picked up.
    void InitScenario();
    /// Record the truth/estimate/covariance for epoch `k` into `res_`.
    void LogEpoch(int k);

    SurfaceRoverNavAppParams params_;
    double t_ = 0.0;

    // Nominal state.
    Vec3d r_ = Vec3d::Zero();
    Vec3d v_ = Vec3d::Zero();
    Mat3d R_b2n_ = Mat3d::Identity();
    Vec3d ba_ = Vec3d::Zero();
    Vec3d bg_ = Vec3d::Zero();
    double cb_ = 0.0;
    double cd_ = 0.0;

    MatXd P_;  // error-state covariance

    // ---- Self-driving (config-constructed) scenario state ---------------------
    bool self_driving_ = false;
    bool initialized_ = false;
    SurfaceNavConfig cfg_;
    SurfaceNavResults res_;

    // RNG shared across Setup and Step so the draw order matches the legacy monolith exactly
    // (a single mt19937 + single normal_distribution, whose Box-Muller cache must persist).
    std::mt19937 rng_;
    std::normal_distribution<double> nd_{0.0, 1.0};
    double Gauss(double sigma) { return sigma * nd_(rng_); }
    Vec3d Gauss3(double sigma) { return Vec3d(Gauss(sigma), Gauss(sigma), Gauss(sigma)); }

    // Local ENU tangent frame (from the World), Moon-fixed geometry.
    int N_ = 0;
    int n_sat_ = 0;
    double dt_ = 1.0;
    double ds_ = 1.0;  // DEM finite-difference step [m]
    Mat3d R_enu2pa_ = Mat3d::Identity();
    Mat3d R_pa2enu_ = Mat3d::Identity();
    Vec3d r_center_pa_ = Vec3d::Zero();
    Vec3d up_hat_pa_ = Vec3d::UnitZ();

    // Precomputed truth series (length N_).
    std::vector<double> Ee_, Nn_, Uu_;
    std::vector<Vec3d> r_truth_, v_truth_;
    std::vector<Mat3d> R_truth_;
    std::vector<Vec3d> f_body_truth_, w_body_truth_;
    std::vector<Vec3d> ba_truth_k_, bg_truth_k_;
    std::vector<std::vector<Vec3d>> sat_pa_;  // [k][j] relay position, Moon-fixed [m]
    std::vector<double> sise_bias_;

    double ClockBiasTruth(int k) const {
      return cfg_.rover_clock_bias_s + cfg_.rover_clock_drift_sps * (k * dt_);
    }
  };

}  // namespace lupnt
