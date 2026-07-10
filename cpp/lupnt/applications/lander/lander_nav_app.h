#pragma once

#include <Eigen/Geometry>
#include <random>
#include <string>
#include <vector>

#include "lupnt/applications/application.h"
#include "lupnt/core/config.h"
#include "lupnt/core/definitions.h"
#include "lupnt/measurements/lander_measurements.h"
#include "lupnt/measurements/surface_measurements.h"

namespace lupnt {

  class LanderGncApp;  // co-hosted guidance app that owns the descent truth trajectory

  /// @brief Onboard tuning for `LanderNavApp` (Kalibr IMU noise densities + measurement /
  /// clock process noise). Mirrors the surface-rover tuning but adds the altimeter and
  /// crater-bearing measurement noise.
  struct LanderNavAppParams {
    // ---- IMU (Kalibr) noise densities ----------------------------------------
    double accel_noise_density = 3.4e-4;  ///< accel white noise [m/s^2 / sqrt(Hz)]
    double accel_bias_rw = 1.0e-4;        ///< accel bias random walk [m/s^3 / sqrt(Hz)]
    double gyro_noise_density = 3.4e-6;   ///< gyro white noise [rad/s / sqrt(Hz)]
    double gyro_bias_rw = 1.0e-6;         ///< gyro bias random walk [rad/s^2 / sqrt(Hz)]

    // ---- Aiding measurement noise --------------------------------------------
    double pseudorange_sigma_m = 1.0;  ///< LunaNet receiver pseudorange noise 1-sigma [m]
    double sise_m = 3.0;               ///< LunaNet signal-in-space error 1-sigma [m]
    double altimeter_sigma_m = 2.0;    ///< radar-altimeter noise 1-sigma [m]
    double crater_sigma_rad = 1.0e-3;  ///< crater-bearing noise 1-sigma [rad]

    // ---- Clock process noise --------------------------------------------------
    double clock_bias_process_sigma = 1.0e-11;   ///< clock-bias random walk [s / sqrt(s)]
    double clock_drift_process_sigma = 1.0e-13;  ///< clock-drift random walk [s/s / sqrt(s)]
  };

  /// @brief Scenario configuration for the lunar-lander navigation application: a lander on
  /// powered descent to a south-pole site, fusing a full IMU (Kalibr), a nadir radar altimeter,
  /// crater-landmark bearings (terrain-relative nav), and LunaNet (LANS) pseudoranges, with the
  /// IMU biases estimated online. The site / DEM is owned by the shared `World`.
  struct LanderNavConfig {
    int seed = 42;
    std::string start_epoch_utc = "2027-03-01T00:00:00";
    double dt_s = 0.5;  ///< IMU / filter / measurement step [s] (must match the guidance app).

    // Site / DEM (site/half_width used only by the legacy free function; the agent-based app
    // takes the DEM from the shared `World`. `dem_max_res_m` is also the terrain-slope FD step).
    double site_lat_deg = -89.45;
    double site_lon_deg = 222.8;
    double dem_half_width_m = 4000.0;
    double dem_max_res_m = 20.0;

    // ---- Lander clock truth (estimated online by the filter) ------------------
    // The descent truth trajectory itself is owned by the co-hosted `LanderGncApp`.
    double lander_clock_bias_s = 1.0e-6;
    double lander_clock_drift_sps = 1.0e-11;

    // ---- IMU (Kalibr noise model) --------------------------------------------
    double accel_noise_density = 3.4e-4;
    double accel_bias_rw = 1.0e-4;
    double gyro_noise_density = 3.4e-6;
    double gyro_bias_rw = 1.0e-6;
    double accel_bias0 = 1.0e-2;
    double gyro_bias0 = 5.0e-4;

    // ---- Radar altimeter ------------------------------------------------------
    bool enable_altimeter = true;
    double altimeter_sigma_m = 2.0;
    double altimeter_max_range_m = 3000.0;

    // ---- Crater-landmark camera (terrain-relative navigation) -----------------
    bool enable_craters = true;
    int n_craters = 60;
    double crater_field_radius_m = 3000.0;
    double camera_fov_deg = 30.0;
    int max_craters_per_epoch = 6;
    double crater_sigma_arcsec = 120.0;

    // ---- LunaNet (LCRNS) constellation ----------------------------------------
    bool enable_lunanet = true;
    std::vector<LcrnsSatConfig> satellites;
    double elevation_mask_deg = 5.0;
    double pseudorange_sigma_m = 1.0;
    double sise_m = 3.0;

    // ---- Filter initialization & tuning --------------------------------------
    double init_pos_sigma_m = 100.0;
    double init_vel_sigma_mps = 1.0;
    double init_att_sigma_deg = 2.0;
    double init_accel_bias_sigma = 2.0e-2;
    double init_gyro_bias_sigma = 1.0e-3;
    double init_clock_bias_sigma_s = 1.0e-6;
    double init_clock_drift_sigma_sps = 1.0e-9;
    double filter_bias_rw_scale = 3.0;
  };

  /// @brief Per-epoch truth/estimate error and covariance series produced by the lander
  /// navigation application, plus the terrain / crater map, in a form directly plottable from
  /// Python. Mirrors the fields the old `RunLanderNav` returned.
  struct LanderNavResults {
    std::string site_id;
    std::string site_name;
    double site_lat_deg = 0.0;
    double site_lon_deg = 0.0;

    MatXd dem_x;
    MatXd dem_y;
    MatXd dem_elevation;
    double dem_center_x = 0.0;
    double dem_center_y = 0.0;
    MatXd crater_enu;  ///< M x 2, crater (East, North) [m].

    VecXd time_s;
    MatXd pos_err_enu;       ///< N x 3, (truth - est) position error in ENU [m].
    MatXd pos_sigma_enu;     ///< N x 3, 1-sigma position uncertainty in ENU [m].
    VecXd pos_err_norm;      ///< N, 3D position error magnitude [m].
    VecXd vel_err_norm;      ///< N, 3D velocity error magnitude [m/s].
    VecXd clock_bias_err;    ///< N, clock-bias error [s].
    VecXd clock_bias_sigma;  ///< N, clock-bias 1-sigma [s].
    VecXi n_visible_sat;     ///< N, number of visible LunaNet satellites.
    VecXi n_craters;         ///< N, number of tracked crater landmarks.

    MatXd accel_bias_err;
    MatXd accel_bias_sigma;
    MatXd gyro_bias_err;
    MatXd gyro_bias_sigma;
    MatXd att_err_deg;
    MatXd att_sigma_deg;
    MatXd att_quat_est;  ///< N x 4, estimated body-to-nav attitude quaternion [w,x,y,z] (MEKF).

    MatXd traj_enu_truth;  ///< N x 3, truth (East, North, Up) [m].
    MatXd traj_enu_est;    ///< N x 3, estimated (East, North, Up) [m].
    VecXd alt_truth;       ///< N, truth height above terrain [m].
    VecXd alt_est;         ///< N, estimated height above terrain [m].

    std::vector<std::string> satellite_names;
  };

  /// @brief Strapdown inertial-navigation **multiplicative EKF (MEKF)** for a lunar lander,
  /// **hosted on a `Lander` agent as its `Application`**. Attitude is estimated as a unit
  /// quaternion `q_b2n` with a 3-parameter multiplicative error in the covariance.
  ///
  /// When constructed from a YAML `application:` block, the app is self-driving: `Setup()`
  /// precomputes the descent truth trajectory (built-in smoothstep or a supplied reference
  /// trajectory), a synthetic crater map, the relay orbits, and the perturbed initial estimate,
  /// seeds the filter, and schedules `Step(t)` at the IMU cadence; each `Step` synthesizes that
  /// epoch's IMU, pseudorange, altimeter, and crater-bearing measurements from the shared
  /// `World` and runs the predict/update cycle, recording the result series accessible below.
  /// The lower-level `Configure`/`Predict`/`Update*` core can also be driven directly.
  class LanderNavApp : public Application {
  public:
    LanderNavApp() = default;
    explicit LanderNavApp(const LanderNavAppParams& params) : params_(params) {}

    /// @brief Construct a self-driving app from the `application:` block of a `Lander` agent.
    explicit LanderNavApp(Config& config);

    /// @brief Precompute descent truth, crater map, relay orbits, and the perturbed initial
    /// estimate (same RNG draw order as the legacy `RunLanderNav`), seed the filter, log epoch 0,
    /// and schedule periodic `Step`s. Terrain / gravity / ENU come from `GetWorld()`.
    void Setup() override;

    /// @brief One self-driving descent epoch: set the host lander's truth state, `Predict` the
    /// IMU sample, then LunaNet / altimeter / crater updates (same order as the legacy loop).
    void Step(Real t) override;

    void Log(Real t) override;

    /// @brief Seed the filter nominal state and covariance.
    void Configure(double t0, const Vec3d& r0, const Vec3d& v0, const Mat3d& R0, const Vec3d& ba0,
                   const Vec3d& bg0, double cb0, double cd0, const MatXd& P0);

    void Predict(const SurfaceImuMeasurement& imu, double dt);
    void UpdateLans(const std::vector<SurfaceLansMeasurement>& meas);
    void UpdateAltimeter(const LanderAltimeterMeasurement& meas);
    void UpdateCraters(const std::vector<LanderCraterMeasurement>& meas);
    void UpdateCrater(const LanderCraterMeasurement& meas);
    void UpdateScalar(const VecXd& H, double z_pred, double z_meas, double variance);

    double time() const { return t_; }
    const Vec3d& position() const { return r_; }
    const Vec3d& velocity() const { return v_; }
    Mat3d attitude() const { return q_b2n_.toRotationMatrix(); }
    Vec4d quaternion() const { return Vec4d(q_b2n_.w(), q_b2n_.x(), q_b2n_.y(), q_b2n_.z()); }
    const Vec3d& accel_bias() const { return ba_; }
    const Vec3d& gyro_bias() const { return bg_; }
    double clock_bias() const { return cb_; }
    double clock_drift() const { return cd_; }
    const MatXd& covariance() const { return P_; }
    const LanderNavAppParams& params() const { return params_; }

    // ---- Result series (valid after the Simulation has run) -------------------
    const LanderNavConfig& config() const { return cfg_; }
    const LanderNavResults& results() const { return res_; }
    const std::string& site_id() const { return res_.site_id; }
    const std::string& site_name() const { return res_.site_name; }
    const MatXd& dem_x() const { return res_.dem_x; }
    const MatXd& dem_y() const { return res_.dem_y; }
    const MatXd& dem_elevation() const { return res_.dem_elevation; }
    const MatXd& crater_enu() const { return res_.crater_enu; }
    const VecXd& time_series() const { return res_.time_s; }
    const MatXd& pos_err_enu() const { return res_.pos_err_enu; }
    const MatXd& pos_sigma_enu() const { return res_.pos_sigma_enu; }
    const VecXd& pos_err_norm() const { return res_.pos_err_norm; }
    const VecXd& vel_err_norm() const { return res_.vel_err_norm; }
    const VecXd& clock_bias_err() const { return res_.clock_bias_err; }
    const VecXd& clock_bias_sigma() const { return res_.clock_bias_sigma; }
    const VecXi& n_visible_sat() const { return res_.n_visible_sat; }
    const VecXi& n_craters() const { return res_.n_craters; }
    const MatXd& accel_bias_err() const { return res_.accel_bias_err; }
    const MatXd& accel_bias_sigma() const { return res_.accel_bias_sigma; }
    const MatXd& gyro_bias_err() const { return res_.gyro_bias_err; }
    const MatXd& gyro_bias_sigma() const { return res_.gyro_bias_sigma; }
    const MatXd& att_err_deg() const { return res_.att_err_deg; }
    const MatXd& att_sigma_deg() const { return res_.att_sigma_deg; }
    const MatXd& att_quat_est() const { return res_.att_quat_est; }
    const MatXd& traj_enu_truth() const { return res_.traj_enu_truth; }
    const MatXd& traj_enu_est() const { return res_.traj_enu_est; }
    const VecXd& alt_truth() const { return res_.alt_truth; }
    const VecXd& alt_est() const { return res_.alt_est; }
    const std::vector<std::string>& satellite_names() const { return res_.satellite_names; }

  private:
    MatXd ProcessNoise(double dt) const;
    void InjectErrorState(const VecXd& dx);
    void UpdateVector(const MatXd& H, const VecXd& y, const MatXd& R);
    /// Precompute descent truth, crater map, relay orbits, and the perturbed initial estimate;
    /// seed the filter and log epoch 0. Deferred to the first `Step` so a reference trajectory
    /// set via `SetReferenceTrajectoryEnu` before `run()` is picked up.
    void InitScenario();
    void LogEpoch(int k);

    LanderNavAppParams params_;
    double t_ = 0.0;

    // Nominal state (attitude as a unit quaternion — the MEKF global orientation).
    Vec3d r_ = Vec3d::Zero();
    Vec3d v_ = Vec3d::Zero();
    Eigen::Quaterniond q_b2n_ = Eigen::Quaterniond::Identity();
    Vec3d ba_ = Vec3d::Zero();
    Vec3d bg_ = Vec3d::Zero();
    double cb_ = 0.0;
    double cd_ = 0.0;

    MatXd P_;  // error-state covariance

    // ---- Self-driving (config-constructed) scenario state ---------------------
    bool self_driving_ = false;
    bool initialized_ = false;
    LanderNavConfig cfg_;
    LanderNavResults res_;

    std::mt19937 rng_;
    std::normal_distribution<double> nd_{0.0, 1.0};
    std::uniform_real_distribution<double> ud_{0.0, 1.0};
    double Gauss(double sigma) { return sigma * nd_(rng_); }
    Vec3d Gauss3(double sigma) { return Vec3d(Gauss(sigma), Gauss(sigma), Gauss(sigma)); }

    int N_ = 0;
    int n_sat_ = 0;
    int n_crat_ = 0;
    double dt_ = 0.5;
    double ds_ = 1.0;
    double crater_sigma_rad_ = 1.0e-3;
    Mat3d R_enu2pa_ = Mat3d::Identity();
    Mat3d R_pa2enu_ = Mat3d::Identity();
    Vec3d r_center_pa_ = Vec3d::Zero();
    Vec3d up_hat_pa_ = Vec3d::UnitZ();

    // The descent truth trajectory is owned by the co-hosted guidance app; read through this
    // handle (resolved in InitScenario). The nav app owns only the sensor/filter state below.
    LanderGncApp* gnc_ = nullptr;
    std::vector<Vec3d> ba_truth_k_, bg_truth_k_;  // truth IMU biases (this app's sensor model)
    std::vector<std::vector<Vec3d>> sat_pa_;
    std::vector<double> sise_bias_;
    std::vector<Vec3d> crater_pa_;

    double ClockBiasTruth(int k) const {
      return cfg_.lander_clock_bias_s + cfg_.lander_clock_drift_sps * (k * dt_);
    }
  };

}  // namespace lupnt
