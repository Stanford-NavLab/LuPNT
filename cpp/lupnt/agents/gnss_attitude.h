/**
 * @file gnss_attitude.h
 * @author Stanford NAV LAB
 * @brief GNSS satellite attitude (body) frame computation
 * @version 0.1
 * @date 2025-06-07
 *
 * @copyright Copyright (c) 2025
 *
 * Computes the canonical GNSS satellite "yaw-steering" body frame used for
 * antenna gain-pattern lookups (nadir / sun cross-track / along-track),
 * mirroring the attitude computation in `GNSSMeas.setup_measurements`
 * (`python/pylupnt/measurements/gnss_meas.py`). This frame is orthonormal by
 * construction and is distinct from the (non-orthonormal) `ijk_to_ecef_rot`
 * frame used for antenna phase-center-offset (PCO) corrections -- see
 * `lupnt/interfaces/antex_loader.h`.
 */
#pragma once

#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief GNSS satellite body (attitude) frame.
  ///
  /// Defines the orthonormal triad `(ex, ey, ez)`:
  ///   - `ez` : nadir direction (toward the central body), `normalize(-r_sat)`
  ///   - `ey` : cross-track direction (toward the Sun side),
  ///            `normalize(cross(ez, normalize(r_sun - r_sat)))`
  ///   - `ex` : along-track direction, completing the right-handed triad,
  ///            `normalize(cross(ey, ez))`
  ///
  /// This is the frame in which transmit-antenna gain patterns are defined
  /// (`theta` = azimuth about `ez` measured from `ex`, `phi` = off-boresight
  /// angle from `ez`).
  class GnssAttitude {
  public:
    GnssAttitude() = default;

    /// @brief Construct and immediately compute the attitude frame.
    ///
    /// Convenience constructor wrapping `Update(r_sat_eci, r_sun_eci)`; used
    /// where a `GnssAttitude` instance (rather than raw `(ex, ey, ez)`
    /// out-parameters) is convenient, e.g. in `test_gnss_attitude.cc`.
    /// @param r_sat_eci Satellite position (ECI or any inertial frame) [m]
    /// @param r_sun_eci Sun position, in the same frame as `r_sat_eci` [m]
    GnssAttitude(const Vec3& r_sat_eci, const Vec3& r_sun_eci);

    /// @brief Compute the canonical GNSS satellite body (attitude) triad
    /// `(ex, ey, ez)` directly from the satellite-to-Sun geometry.
    ///
    /// This is the core attitude computation used to evaluate transmit
    /// antenna gain patterns: `GNSSMeasurements::ComputeCN0`
    /// (`lupnt/measurements/gnss_measurement.cc`) calls the
    /// velocity-aware overload below to get the transmitter's `(ex, ey, ez)`
    /// at the signal-transmission epoch, then projects the
    /// transmitter-to-receiver line of sight onto this triad to get the
    /// off-boresight angles fed to `Antenna::ComputeGain`. Mirrors the
    /// attitude computation in `GNSSMeas.setup_measurements`
    /// (`python/pylupnt/measurements/gnss_meas.py`).
    /// @param r_sat_eci Satellite position (ECI or any inertial frame) [m]
    /// @param r_sun_eci Sun position, in the same frame as `r_sat_eci` [m]
    /// @param[out] ex Along-track body axis (unit vector)
    /// @param[out] ey Cross-track body axis (unit vector)
    /// @param[out] ez Nadir body axis (unit vector)
    static void Compute(const Vec3& r_sat_eci, const Vec3& r_sun_eci, Vec3& ex, Vec3& ey, Vec3& ez);

    /// @brief Compute the attitude triad `(ex, ey, ez)` via the documented
    /// nominal yaw-steering law (`GnssYawSteering::NominalYawAngle`, Eq. 1 of
    /// Cheng et al., 2025), rather than directly from the Sun direction.
    ///
    /// Internally computes the Sun-elevation (`beta`) and orbit (`mu`) angles
    /// via `GnssYawSteering::BetaAngle`/`OrbitAngle`, evaluates the nominal
    /// yaw angle `phi = NominalYawAngle(beta, mu)`, and builds the frame by
    /// rotating the orbital reference frame about the nadir axis by `phi`
    /// (see `ComputeFromYawAngle`). This is *numerically identical* to
    /// `Compute(r_sat_eci, r_sun_eci, ...)` -- both describe the same
    /// Sun-pointing nominal attitude, since the nominal yaw law is, by
    /// construction, the angle that keeps the solar panels Sun-pointing --
    /// but this overload makes the dependency on the documented yaw-steering
    /// law explicit, and `ComputeFromYawAngle` is the basis for computing the
    /// *modeled* (maneuver/eclipse) attitude from any of `GnssYawSteering`'s
    /// block-specific yaw-steering laws (Eq. 3-7, 12-16). This is the overload
    /// called by `GNSSMeasurements::ComputeCN0` to get the transmitter's body
    /// frame for antenna-gain lookups.
    /// @param r_sat_eci Satellite position (ECI or any inertial frame) [m]
    /// @param v_sat_eci Satellite velocity, in the same frame [m/s]
    /// @param r_sun_eci Sun position, in the same frame as `r_sat_eci` [m]
    /// @param[out] ex Along-track-ish body axis (unit vector)
    /// @param[out] ey Cross-track-ish body axis (unit vector)
    /// @param[out] ez Nadir body axis (unit vector)
    static void Compute(const Vec3& r_sat_eci, const Vec3& v_sat_eci, const Vec3& r_sun_eci,
                        Vec3& ex, Vec3& ey, Vec3& ez);

    /// @brief Build the attitude triad `(ex, ey, ez)` from an explicit yaw
    /// angle `phi`, e.g. one produced by any of `GnssYawSteering`'s nominal
    /// or modeled/maneuver yaw-steering laws (Eq. 1, 3-7, 8, 12-16), and the
    /// orbital geometry alone (no Sun direction needed):
    ///   `ez     = normalize(-r_sat)`                       (nadir)
    ///   `ex_ref = normalize(v_sat)`, `ey_ref = ez x ex_ref`   (orbital frame, phi = 0)
    ///   `ex = cos(phi) * ex_ref + sin(phi) * ey_ref`
    ///   `ey = -sin(phi) * ex_ref + cos(phi) * ey_ref`
    /// i.e. `(ex, ey)` are the orbital reference axes `(ex_ref, ey_ref)`
    /// rotated by `phi` about the nadir axis `ez`. Supplying
    /// `phi = GnssYawSteering::NominalYawAngle(beta, mu)` reproduces exactly
    /// the Sun-pointing nominal frame of `Compute` (verified in
    /// `test_gnss_attitude.cc`); supplying a *modeled* yaw angle from one of
    /// `GnssYawSteering`'s maneuver-specific laws instead yields the actual
    /// (eclipse-season / midnight-noon-turn) attitude.
    /// @param r_sat_eci Satellite position (ECI or any inertial frame) [m]
    /// @param v_sat_eci Satellite velocity, in the same frame [m/s]
    /// @param yaw_angle Yaw angle `phi` [rad]
    /// @param[out] ex Along-track-ish body axis (unit vector)
    /// @param[out] ey Cross-track-ish body axis (unit vector)
    /// @param[out] ez Nadir body axis (unit vector)
    static void ComputeFromYawAngle(const Vec3& r_sat_eci, const Vec3& v_sat_eci, Real yaw_angle,
                                    Vec3& ex, Vec3& ey, Vec3& ez);

    /// @brief (Re-)compute and cache the attitude triad for this instance.
    void Update(const Vec3& r_sat_eci, const Vec3& r_sun_eci);

    /// @brief (Re-)compute and cache the attitude triad via the documented
    /// nominal yaw-steering law; see
    /// `Compute(r_sat_eci, v_sat_eci, r_sun_eci, ...)`.
    void Update(const Vec3& r_sat_eci, const Vec3& v_sat_eci, const Vec3& r_sun_eci);

    /// @brief Get the cached along-track body axis `ex` (unit vector, in the
    /// same inertial frame as the `r_sat_eci`/`r_sun_eci` passed to
    /// `Update`/`Compute`). Part of this orthonormal, GNSS-attitude-frame
    /// triad -- distinct from the non-orthonormal `ijk` PCO frame of
    /// `AntexLoader::ComputeIjkToEcefRotation`.
    const Vec3& GetEx() const { return ex_; }
    /// @brief Get the cached cross-track (Sun-side) body axis `ey` (unit
    /// vector, same frame as `GetEx`). Part of this orthonormal,
    /// GNSS-attitude-frame triad -- distinct from the non-orthonormal `ijk`
    /// PCO frame of `AntexLoader::ComputeIjkToEcefRotation`.
    const Vec3& GetEy() const { return ey_; }
    /// @brief Get the cached nadir body axis `ez` (unit vector, same frame as
    /// `GetEx`). Part of this orthonormal, GNSS-attitude-frame triad --
    /// distinct from the non-orthonormal `ijk` PCO frame of
    /// `AntexLoader::ComputeIjkToEcefRotation`.
    const Vec3& GetEz() const { return ez_; }

    /// @brief Body-to-inertial rotation matrix `[ex, ey, ez]` (columns),
    /// assembled from the cached triad. Used to express vectors given in this
    /// GNSS-attitude body frame (e.g. antenna boresight/PCO offsets) in the
    /// same inertial frame as `r_sat_eci`.
    Mat3 GetRotationMatrix() const;

    /// @brief Off-boresight angle `theta` (azimuth, about `ez` measured from
    /// `ex`) and `phi` (polar angle from `ez`) of a unit direction `u`
    /// (expressed in the same inertial frame as `r_sat_eci`/`r_sun_eci`),
    /// as seen from the satellite body frame.
    ///
    /// Mirrors the angle computation used by GNSS link-budget and
    /// measurement-channel construction: given the line-of-sight direction
    /// to a receiver/target, returns the `(theta, phi)` pair that
    /// `Antenna::ComputeGain` expects to look up the transmit/receive gain
    /// pattern.
    /// @param u     Unit direction vector, in the same inertial frame as
    ///        `r_sat_eci`/`r_sun_eci`
    /// @param[out] theta Azimuth angle about `ez`, measured from `ex` [rad]
    /// @param[out] phi   Polar (off-boresight) angle from `ez` [rad]
    void GetAngles(const Vec3& u, Real& theta, Real& phi) const;

  private:
    Vec3 ex_ = Vec3::UnitX();
    Vec3 ey_ = Vec3::UnitY();
    Vec3 ez_ = Vec3::UnitZ();
  };

}  // namespace lupnt
