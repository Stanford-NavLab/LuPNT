#pragma once

#include <string>

#include "lupnt/core/definitions.h"
#include "lupnt/measurements/measurement.h"
#include "lupnt/measurements/measurement_utils.h"

namespace lupnt {

  /// @brief Initial Cartesian state of one LCRNS / LunaNet relay satellite, defined at the
  /// scenario start epoch in the Moon-centered inertial frame (`Frame::MOON_CI`), e.g. taken
  /// from the NASA LCRNS reference constellation. Shared by the surface-rover and lander
  /// navigation applications (relays are propagated as Keplerian orbits app-internally).
  struct LcrnsSatConfig {
    std::string name = "SV";
    Vec3d r0_m = Vec3d::Zero();    ///< initial position [m], Frame::MOON_CI.
    Vec3d v0_mps = Vec3d::Zero();  ///< initial velocity [m/s], Frame::MOON_CI.
  };

  /// Error-state dimension of the surface-navigation INS filter:
  /// `[dr(3), dv(3), dtheta(3), db_a(3), db_g(3), d(clock_bias), d(clock_drift)]`.
  constexpr int kSurfaceNavErrorStateSize = 17;

  /// @brief One inertial-measurement-unit sample used by the strapdown INS.
  ///
  /// Carries the raw **body-frame** accelerometer specific force (`accel`) and gyroscope
  /// angular rate (`gyro`) -- i.e. the true quantities corrupted by a slowly-varying bias
  /// and white noise, following the Kalibr IMU noise model
  /// (https://github.com/ethz-asl/kalibr/wiki/IMU-Noise-Model). The filter corrects these
  /// with its estimated biases and mechanizes them into position/velocity/attitude, so no
  /// helper methods live here -- the mechanization needs the full nav state and attitude
  /// (see `SurfaceRoverNavApp::Predict`).
  struct SurfaceImuMeasurement {
    double timestamp = 0.0;       ///< sample epoch [s]
    Vec3d accel = Vec3d::Zero();  ///< body-frame specific force [m/s^2]
    Vec3d gyro = Vec3d::Zero();   ///< body-frame angular rate [rad/s]
  };

  /// @brief One LCRNS / LANS (Lunar Augmented Navigation Service) pseudorange from a relay
  /// satellite to the surface rover.
  ///
  /// Models a one-way code pseudorange biased by the rover clock:
  /// `rho = |r_sat - r_rover| + C * clock_bias`. The measurement noise combines the
  /// receiver pseudorange noise (`sigma_m`, default 1 m) and the broadcast ephemeris/clock
  /// signal-in-space error (`sise_m`, default 3 m, 1-sigma) added in quadrature.
  struct SurfaceLansMeasurement : public ErrorStateMeasurement {
    /// @brief Error-state mapping for the LANS pseudorange model.
    struct Config {
      int i_dr = 0;    ///< position-error offset in the error state
      int i_dcb = 15;  ///< clock-bias-error offset in the error state
    };

    double timestamp = 0.0;       ///< measurement epoch [s]
    Vec3d r_sat = Vec3d::Zero();  ///< transmitter position, Moon-fixed frame [m]
    double pseudorange_m = 0.0;   ///< measured pseudorange [m]
    double sigma_m = 1.0;         ///< receiver pseudorange noise 1-sigma [m]
    double sise_m = 3.0;          ///< signal-in-space error 1-sigma [m]
    Config config;                ///< error-state index mapping

    /// @brief Predicted pseudorange `|r_sat - r_rover| + C * clock_bias`.
    /// @param r_rover     Rover position, Moon-fixed frame [m].
    /// @param clock_bias_s Rover clock bias [s].
    double PredictedRange(const Vec3d& r_rover, double clock_bias_s) const;

    /// @brief Unit line-of-sight vector from the rover toward the satellite.
    /// The range partial `d(rho)/d(r_rover)` equals `-LosUnit`.
    /// @param r_rover Rover position, Moon-fixed frame [m].
    Vec3d LosUnit(const Vec3d& r_rover) const;

    /// @brief Measurement noise variance `sigma_m^2 + sise_m^2` [m^2].
    double NoiseVariance() const { return sigma_m * sigma_m + sise_m * sise_m; }

    /// @brief Error-state model: predicted pseudorange, `R`, and (if `H != nullptr`) the
    /// `1 x error_state_size` Jacobian (`-LosUnit^T` at the position error, `C` at the
    /// clock-bias error). Consumes only the nominal position and clock bias.
    MeasData Compute(const NavErrorContext& nom, MatXd* H = nullptr) const override;
  };

}  // namespace lupnt
