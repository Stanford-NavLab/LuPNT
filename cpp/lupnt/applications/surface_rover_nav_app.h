#pragma once

#include <vector>

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

  /// @brief Strapdown inertial-navigation error-state EKF that fuses a full IMU
  /// (accelerometer + gyroscope), LCRNS (LANS) pseudoranges, and a DEM altitude constraint
  /// to position and orient a surface rover, **estimating the IMU biases online**.
  ///
  /// The nominal state is `[r(3), v(3), R_b2n(SO3), b_a(3), b_g(3), clock_bias, clock_drift]`
  /// in a Moon-fixed frame; the 17-element error state
  /// `[dr, dv, dtheta, db_a, db_g, d(cb), d(cd)]` (see `kSurfaceNavErrorStateSize`) carries
  /// the covariance. Each step the driver `Predict`s with an IMU sample (bias-corrected
  /// strapdown mechanization + Kalibr process noise), then applies scalar updates: one per
  /// visible LANS satellite (`UpdateLans`) and, optionally, a linearized DEM altitude
  /// constraint (`UpdateScalar`). Attitude observability comes from the accelerometer
  /// sensing the gravity direction plus rover motion; gyro-bias observability follows from
  /// the attitude, which is why the full attitude state is carried.
  ///
  /// Kept independent of the `Application`/`Agent` scheduler so it can be driven directly by
  /// `SurfaceNavSimulation` (and unit-tested).
  class SurfaceRoverNavApp {
  public:
    SurfaceRoverNavApp() = default;
    explicit SurfaceRoverNavApp(const SurfaceRoverNavAppParams& params) : params_(params) {}

    /// @brief Seed the filter nominal state and covariance.
    /// @param t0    Initial epoch [s].
    /// @param r0    Initial position [m], Moon-fixed frame.
    /// @param v0    Initial velocity [m/s], Moon-fixed frame.
    /// @param R0    Initial body-to-nav rotation (SO3).
    /// @param ba0   Initial accelerometer bias estimate [m/s^2], body frame.
    /// @param bg0   Initial gyroscope bias estimate [rad/s], body frame.
    /// @param cb0   Initial clock bias [s].
    /// @param cd0   Initial clock drift [s/s].
    /// @param P0    Initial error-state covariance (`kSurfaceNavErrorStateSize` square).
    void Configure(double t0, const Vec3d& r0, const Vec3d& v0, const Mat3d& R0, const Vec3d& ba0,
                   const Vec3d& bg0, double cb0, double cd0, const MatXd& P0);

    /// @brief Predict to `t = t_ + dt` using a full IMU sample: bias-correct the raw
    /// accelerometer/gyro, mechanize `R_b2n`, `v`, `r` (adding lunar gravity), advance the
    /// clock, and propagate the error-state covariance with the INS transition matrix and
    /// Kalibr process noise.
    /// @param imu IMU sample (body-frame specific force + angular rate).
    /// @param dt  Time step [s].
    void Predict(const SurfaceImuMeasurement& imu, double dt);

    /// @brief Apply a sequential scalar EKF update for each LANS pseudorange.
    void UpdateLans(const std::vector<SurfaceLansMeasurement>& meas);

    /// @brief Generic scalar EKF update in the error state, immediately injecting the
    /// resulting correction into the nominal state (error-state reset).
    /// @param H        Measurement Jacobian wrt the error state (`kSurfaceNavErrorStateSize`).
    /// @param z_pred   Predicted measurement.
    /// @param z_meas   Observed measurement.
    /// @param variance Measurement noise variance.
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

  private:
    /// Build the discrete process-noise covariance for a step of length `dt`.
    MatXd ProcessNoise(double dt) const;
    /// Inject a full error-state correction into the nominal state.
    void InjectErrorState(const VecXd& dx);

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
  };

}  // namespace lupnt
