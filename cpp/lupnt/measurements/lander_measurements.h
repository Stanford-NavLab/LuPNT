#pragma once

#include <string>

#include "lupnt/core/definitions.h"
#include "lupnt/measurements/measurement.h"
#include "lupnt/measurements/measurement_utils.h"
#include "lupnt/measurements/surface_measurements.h"

namespace lupnt {

  /// Error-state dimension of the lunar-lander INS filter:
  /// `[dr(3), dv(3), dtheta(3), db_a(3), db_g(3), d(clock_bias), d(clock_drift)]`.
  /// Identical layout to the surface-rover filter (`kSurfaceNavErrorStateSize`).
  constexpr int kLanderNavErrorStateSize = 17;

  // The lander reuses the surface IMU and LANS/LunaNet pseudorange measurement models:
  //   - `SurfaceImuMeasurement` : body-frame accelerometer specific force + gyro angular rate.
  //   - `SurfaceLansMeasurement`: one-way LunaNet (LANS) pseudorange from a relay satellite.
  // (see `lupnt/measurements/surface_measurements.h`).

  /// @brief One radar/laser altimeter return: the vehicle's height above the terrain
  /// directly below it (nadir), i.e. `altitude = U_vehicle - terrain(East, North)`.
  ///
  /// Modelled as a vertical altitude measurement (small-tilt / nadir-beam assumption). The
  /// measurement Jacobian couples to position through the local East-North-Up rotation and the
  /// terrain slope `(dterrain/dE, dterrain/dN)`; because those depend on the DEM they are
  /// supplied by the simulation, which drives the update through `LanderNavApp::UpdateAltimeter`.
  struct LanderAltimeterMeasurement : public ErrorStateMeasurement {
    /// @brief Error-state mapping for the altimeter model.
    struct Config {
      int i_dr = 0;  ///< position-error offset in the error state
    };

    double timestamp = 0.0;   ///< measurement epoch [s].
    double altitude_m = 0.0;  ///< measured height above terrain [m].
    double sigma_m = 2.0;     ///< altimeter noise 1-sigma [m].

    /// Predicted altitude at the nominal position and the position partial `d(altitude)/d(r)`
    /// are DEM-dependent, so they are supplied by the simulation before the update.
    double predicted_altitude_m = 0.0;  ///< DEM-predicted altitude at nominal position [m].
    Vec3d h_pos = Vec3d::Zero();        ///< d(altitude)/d(r), Moon-fixed frame (DEM slope).
    Config config;                      ///< error-state index mapping.

    /// @brief Measurement noise variance `sigma_m^2` [m^2].
    double NoiseVariance() const { return sigma_m * sigma_m; }

    /// @brief Error-state model: predicted altitude, `R`, and (if `H != nullptr`) the
    /// `1 x error_state_size` Jacobian (`h_pos^T` at the position error).
    MeasData Compute(const NavErrorContext& nom, MatXd* H = nullptr) const override;
  };

  /// @brief One crater-landmark bearing (terrain-relative navigation): the unit line-of-sight
  /// from the lander's down-looking camera to a mapped crater whose Moon-fixed position is
  /// known, expressed in the vehicle **body** frame.
  ///
  /// Predicted body line-of-sight is `R_b2n^T (r_crater - r) / |r_crater - r|`, so the
  /// residual constrains both horizontal position and attitude. The angular noise `sigma_rad`
  /// is applied per line-of-sight component (small-angle unit-vector approximation).
  struct LanderCraterMeasurement : public ErrorStateMeasurement {
    /// @brief Error-state mapping for the crater-bearing model.
    struct Config {
      int i_dr = 0;   ///< position-error offset in the error state
      int i_dth = 6;  ///< attitude-error offset in the error state
    };

    double timestamp = 0.0;          ///< measurement epoch [s].
    Vec3d r_crater = Vec3d::Zero();  ///< known landmark position, Moon-fixed frame [m].
    Vec3d los_body = Vec3d::Zero();  ///< measured unit line-of-sight in the body frame.
    double sigma_rad = 1.0e-3;       ///< bearing noise 1-sigma [rad].
    std::string id;                  ///< landmark identifier.
    Config config;                   ///< error-state index mapping.

    /// @brief Predicted unit line-of-sight in the body frame:
    /// `R_b2n^T (r_crater - r) / |r_crater - r|`.
    /// @param r     Lander position, Moon-fixed frame [m].
    /// @param R_b2n Body-to-nav rotation.
    Vec3d PredictedLosBody(const Vec3d& r, const Mat3d& R_b2n) const;

    /// @brief Error-state model: predicted body line-of-sight, `R = sigma_rad^2 I`, and
    /// (if `H != nullptr`) the `3 x error_state_size` Jacobian coupling position
    /// (`R_b2n^T d(u_n)/d(r)`) and attitude (`R_b2n^T [u_n]_x`). Returns an empty value
    /// vector for a degenerate (zero-range) geometry.
    MeasData Compute(const NavErrorContext& nom, MatXd* H = nullptr) const override;
  };

}  // namespace lupnt
