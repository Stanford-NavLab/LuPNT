#pragma once

#include <Eigen/Geometry>
#include <vector>

#include "lupnt/applications/application.h"
#include "lupnt/core/definitions.h"
#include "lupnt/measurements/lander_measurements.h"
#include "lupnt/measurements/surface_measurements.h"

namespace lupnt {

  /// @brief Onboard tuning for `LanderNavApp` (Kalibr IMU noise densities + measurement /
  /// clock process noise). Mirrors the surface-rover tuning but adds the altimeter and
  /// crater-bearing measurement noise. IMU defaults are representative tactical-grade
  /// (LN200S-like) values (Kalibr model: https://github.com/ethz-asl/kalibr/wiki/IMU-Noise-Model).
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

  /// @brief Strapdown inertial-navigation **multiplicative EKF (MEKF)** for a lunar lander,
  /// **hosted on a `Lander` agent as its `Application`** (attach via `Agent::SetApplication`).
  ///
  /// Attitude is estimated as a **unit quaternion** `q_b2n` (body-to-nav) carried in the nominal
  /// state, with a minimal 3-parameter multiplicative attitude error `dtheta` (nav-frame
  /// rotation vector) in the covariance. This is the defining feature of the MEKF: the global
  /// orientation lives on the quaternion manifold (never over-parameterized in the covariance),
  /// and each update injects its correction *multiplicatively* — `q_b2n <- dq(dtheta) (x) q_b2n`
  /// — then resets the error to zero. Propagation integrates the quaternion kinematics from the
  /// bias-corrected gyro.
  ///
  /// The nominal state is `[r(3), v(3), q_b2n(unit quat), b_a(3), b_g(3), clock_bias,
  /// clock_drift]` in a Moon-fixed frame; the 17-element error state
  /// `[dr, dv, dtheta, db_a, db_g, d(cb), d(cd)]` (see `kLanderNavErrorStateSize`) carries the
  /// covariance. Each descent epoch the simulation `Predict`s with an IMU sample (bias-corrected
  /// strapdown mechanization + Kalibr process noise), then applies aiding updates:
  ///   - `UpdateLans`      : one scalar update per visible LunaNet (LANS) pseudorange;
  ///   - `UpdateAltimeter` : a scalar height-above-terrain update (nadir altimeter);
  ///   - `UpdateCraters`   : a 3-vector body line-of-sight update per mapped crater landmark,
  ///                         constraining horizontal position and attitude (terrain-relative nav).
  ///
  /// The estimation core is deliberately driven directly by `RunLanderNav` (rather than the
  /// `Simulation` event scheduler), so `Step` is a no-op; the value of the `Application` wiring
  /// is that the filter genuinely belongs to — and is reachable from — the lander agent.
  class LanderNavApp : public Application {
  public:
    LanderNavApp() = default;
    explicit LanderNavApp(const LanderNavAppParams& params) : params_(params) {}

    /// @brief No-op: the descent simulation drives `Predict`/`Update*` explicitly rather than
    /// through the periodic scheduler.
    void Step(Real /*t*/) override {}

    /// @brief Seed the filter nominal state and covariance.
    /// @param t0    Initial epoch [s].
    /// @param r0    Initial position [m], Moon-fixed frame.
    /// @param v0    Initial velocity [m/s], Moon-fixed frame.
    /// @param R0    Initial body-to-nav rotation (SO3); stored internally as the unit quaternion.
    /// @param ba0   Initial accelerometer bias estimate [m/s^2], body frame.
    /// @param bg0   Initial gyroscope bias estimate [rad/s], body frame.
    /// @param cb0   Initial clock bias [s].
    /// @param cd0   Initial clock drift [s/s].
    /// @param P0    Initial error-state covariance (`kLanderNavErrorStateSize` square).
    void Configure(double t0, const Vec3d& r0, const Vec3d& v0, const Mat3d& R0, const Vec3d& ba0,
                   const Vec3d& bg0, double cb0, double cd0, const MatXd& P0);

    /// @brief Predict to `t = t_ + dt` using a full IMU sample (bias-correct, mechanize
    /// `R_b2n`, `v`, `r` with lunar gravity, advance the clock, propagate the covariance).
    void Predict(const SurfaceImuMeasurement& imu, double dt);

    /// @brief Apply a sequential scalar EKF update for each LunaNet (LANS) pseudorange.
    void UpdateLans(const std::vector<SurfaceLansMeasurement>& meas);

    /// @brief Apply the nadir radar-altimeter height-above-terrain update.
    /// @param meas Altimeter return; carries the measured height plus the DEM-derived
    ///             `predicted_altitude_m` and position partial `h_pos` (supplied by the
    ///             simulation from the local ENU frame and terrain slope).
    void UpdateAltimeter(const LanderAltimeterMeasurement& meas);

    /// @brief Apply a crater-bearing (body line-of-sight) update per landmark.
    void UpdateCraters(const std::vector<LanderCraterMeasurement>& meas);
    /// @brief Apply one crater-bearing update.
    void UpdateCrater(const LanderCraterMeasurement& meas);

    /// @brief Generic scalar EKF update in the error state, injecting the correction into the
    /// nominal state (error-state reset).
    void UpdateScalar(const VecXd& H, double z_pred, double z_meas, double variance);

    double time() const { return t_; }            ///< current filter epoch [s].
    const Vec3d& position() const { return r_; }  ///< estimated position [m].
    const Vec3d& velocity() const { return v_; }  ///< estimated velocity [m/s].
    /// Estimated body-to-nav rotation (derived from the nominal quaternion).
    Mat3d attitude() const { return q_b2n_.toRotationMatrix(); }
    /// Estimated body-to-nav attitude quaternion as `[w, x, y, z]` (MEKF nominal state).
    Vec4d quaternion() const { return Vec4d(q_b2n_.w(), q_b2n_.x(), q_b2n_.y(), q_b2n_.z()); }
    const Vec3d& accel_bias() const { return ba_; }  ///< estimated accel bias [m/s^2].
    const Vec3d& gyro_bias() const { return bg_; }   ///< estimated gyro bias [rad/s].
    double clock_bias() const { return cb_; }        ///< estimated clock bias [s].
    double clock_drift() const { return cd_; }       ///< estimated clock drift [s/s].
    const MatXd& covariance() const { return P_; }   ///< error-state covariance.
    const LanderNavAppParams& params() const { return params_; }

  private:
    MatXd ProcessNoise(double dt) const;
    void InjectErrorState(const VecXd& dx);
    /// Vector EKF update with measurement Jacobian `H` (m x NX), residual `y` (m), noise `R`.
    void UpdateVector(const MatXd& H, const VecXd& y, const MatXd& R);

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
  };

}  // namespace lupnt
