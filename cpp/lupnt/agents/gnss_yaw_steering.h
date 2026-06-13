/**
 * @file gnss_yaw_steering.h
 * @author Stanford NAV LAB
 * @brief GNSS satellite yaw-attitude steering laws
 * @version 0.1
 * @date 2026-06-07
 *
 * @copyright Copyright (c) 2026
 *
 * Implements the nominal GNSS yaw-attitude model and the precise
 * ("modeled") yaw-steering laws used by GPS, Galileo and BDS-3 satellites
 * during yaw maneuvers (orbit midnight/noon turns and eclipse-season yaw
 * steering), following:
 *
 *   L. Cheng, T. Geng, J. Liu, X. Xie, P. Zhao, W. Liu, "GNSS satellite yaw
 *   attitude laws and impact on satellite clock, phase OSB and PPP-AR",
 *   Advances in Space Research 75 (2025) 2535-2549,
 *   https://doi.org/10.1016/j.asr.2024.10.064
 *
 * Two angles parameterize the satellite-Sun-Earth geometry that drives all
 * of these laws (Kouba, 2009):
 *   - `beta` : Sun elevation angle above the orbital plane [rad] -- positive
 *              when the Sun is on the same side of the orbital plane as the
 *              orbit-normal (angular-momentum) vector.
 *   - `mu`   : orbit angle [rad] -- the geocentric angle, measured in the
 *              orbital plane in the direction of motion, from the orbit
 *              "midnight point" (where the satellite is directly behind the
 *              Earth as seen from the Sun) to the satellite. `mu = 0` at the
 *              midnight point and `mu = +/-pi` at the "noon point".
 *
 * Each law below is a stateless, pure function mirroring one numbered
 * equation from the paper; they are intentionally low-level building blocks
 * (no maneuver-window detection / state machine is implemented, since the
 * trigger conditions differ across satellite blocks and are only described
 * qualitatively in the reference). A typical caller:
 *   1. computes `beta`/`mu` (and their rates) from the satellite ephemeris
 *      and Sun position, e.g. via `BetaAngle`/`OrbitAngle` below,
 *   2. evaluates `NominalYawAngle` to get the nominal attitude,
 *   3. detects whether the satellite is inside a yaw-maneuver window
 *      (shadow crossing / midnight or noon turn / eclipse season) using the
 *      block-specific criteria summarized in the reference (Section 2), and
 *   4. if so, evaluates the corresponding "modeled" law instead.
 */
#pragma once

#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief GNSS satellite yaw-attitude steering laws (Cheng et al., 2025).
  ///
  /// All angles are in radians and all angular rates are in rad/s, in
  /// keeping with the rest of `lupnt` (the empirical constants from the
  /// paper, e.g. `0.11 deg/s`, `3 deg`, are converted internally via `RAD`).
  /// All functions are stateless static members; the class exists purely to
  /// namespace the laws (mirroring `GnssAttitude`/`AntexLoader`).
  class GnssYawSteering {
  public:
    GnssYawSteering() = delete;

    // ---- Satellite-Sun-Earth geometry --------------------------------------

    /// @brief Sun elevation angle `beta` above the orbital plane [rad].
    /// `beta = asin(n . s_hat)`, where `n = (r x v)/|r x v|` is the orbit
    /// normal (angular-momentum direction) and `s_hat = r_sun/|r_sun|`.
    ///
    /// Step 1 of the yaw-steering pipeline described above: called by
    /// `GnssAttitude::Compute(r_sat_eci, v_sat_eci, r_sun_eci, ...)` to derive
    /// the satellite-Sun-Earth geometry needed to evaluate `NominalYawAngle`.
    /// @param r_sat Satellite position (ECI or any inertial frame) [m]
    /// @param v_sat Satellite velocity, in the same frame as `r_sat` [m/s]
    /// @param r_sun Sun position, in the same frame as `r_sat` [m]
    static Real BetaAngle(const Vec3& r_sat, const Vec3& v_sat, const Vec3& r_sun);

    /// @brief Orbit angle `mu` [rad]: the geocentric angle, in the orbital
    /// plane and in the direction of motion, from the orbit "midnight point"
    /// to the satellite (`mu = 0` at the midnight point, `mu = +/-pi` at the
    /// "noon point", `mu` increases in the direction of orbital motion).
    ///
    /// Step 1 of the yaw-steering pipeline described above: called alongside
    /// `BetaAngle` by `GnssAttitude::Compute(r_sat_eci, v_sat_eci, r_sun_eci,
    /// ...)` to derive the satellite-Sun-Earth geometry needed to evaluate
    /// `NominalYawAngle`.
    /// @param r_sat Satellite position (ECI or any inertial frame) [m]
    /// @param v_sat Satellite velocity, in the same frame as `r_sat` [m/s]
    /// @param r_sun Sun position, in the same frame as `r_sat` [m]
    static Real OrbitAngle(const Vec3& r_sat, const Vec3& v_sat, const Vec3& r_sun);

    /// @brief Convert a midnight-point-referenced orbit angle `mu` to the
    /// noon-point-referenced orbit angle `eta` used in Eq. (9): the noon
    /// point is diametrically opposite ( +/- pi ) the midnight point along
    /// the orbit, so `eta = WrapToPi(mu - pi)`.
    ///
    /// Bridges `OrbitAngle`'s midnight-referenced `mu` (used by
    /// `NominalYawAngle`/`GnssAttitude`) to the noon-referenced `eta`
    /// convention used by the Galileo IOV laws below
    /// (`GalileoSunVector`/`GalileoIovNominalYawAngle`/`GalileoIovEclipseYawAngle`).
    /// @param mu Orbit angle from the midnight point [rad]
    /// @return   Orbit angle from the noon point [rad]
    static Real OrbitNoonAngle(Real mu);

    // ---- Nominal yaw-attitude model (Kouba, 2009), Eq. (1)-(2) --------------

    /// @brief Nominal yaw angle                                        [Eq. 1]
    ///   `phi(t) = ATAN2(-tan(beta), sin(mu))`
    /// Eq. (1) shows that the sign of the nominal yaw angle is always
    /// opposite to that of `beta`.
    ///
    /// Step 2 of the yaw-steering pipeline described above: evaluated by
    /// `GnssAttitude::Compute(r_sat_eci, v_sat_eci, r_sun_eci, ...)` from the
    /// `beta`/`mu` computed via `BetaAngle`/`OrbitAngle`, then passed to
    /// `GnssAttitude::ComputeFromYawAngle` to build the Sun-pointing nominal
    /// body frame.
    static Real NominalYawAngle(Real beta, Real mu);

    /// @brief Nominal yaw rate                                         [Eq. 2]
    ///   `phi_dot = mu_dot * tan(beta) * cos(mu) / (sin(mu)^2 + tan(beta)^2)`
    /// where `mu_dot` is the (average) orbital angular rate.
    ///
    /// Provided as a low-level building block for maneuver-window detection
    /// (step 3 of the yaw-steering pipeline described above); not currently
    /// called elsewhere in LuPNT, since no maneuver-window state machine is
    /// implemented yet.
    static Real NominalYawRate(Real beta, Real mu, Real mu_dot);

    // ---- GPS Block IIF -------------------------------------------------------
    // (Kuang et al., 2017)

    /// @brief Constant-yaw-rate model for GPS Block IIF satellites while
    /// crossing the Earth's shadow                                   [Eq. 3-4]
    /// The satellite maintains a constant yaw rate
    ///   `phi_dot_star = (phi(t_e) - phi(t_s)) / (t_e - t_s)`
    /// from the shadow-entry epoch `t_s` to the shadow-exit epoch `t_e`, so
    ///   `phi(t) = phi(t_s) + phi_dot_star * (t - t_s)`.
    /// @param t Evaluation epoch (TAI seconds, or any consistent time base)
    /// @param ts Shadow entry epoch
    /// @param te Shadow exit epoch
    /// @param phi_ts Nominal yaw angle at `ts` [rad]
    /// @param phi_te Nominal yaw angle at `te` [rad]
    static Real GpsIIFShadowYawAngle(Real t, Real ts, Real te, Real phi_ts, Real phi_te);

    /// @brief Noon-turn maneuver model for GPS Block IIF satellites    [Eq. 5]
    ///   `phi(t) = phi(t_s) - SIGN(R, beta + 0.7 deg) * (t - t_s)`
    /// with maximum hardware yaw rate `R = 0.11 deg/s` (Kouba, 2013; Liu
    /// et al., 2018). The `0.7 deg` bias avoids reverse yaw maneuvers when
    /// `beta` falls in `[-0.7, 0] deg` (Kuang et al., 2017).
    /// @param t Evaluation epoch
    /// @param ts Noon-turn start epoch
    /// @param phi_ts Nominal yaw angle at `ts` [rad]
    /// @param beta Sun elevation angle above the orbital plane [rad]
    static Real GpsIIFNoonTurnYawAngle(Real t, Real ts, Real phi_ts, Real beta);

    // ---- GPS Block IIR -------------------------------------------------------
    // (Kouba, 2017)

    /// @brief Midnight-turn maneuver model for GPS Block IIR satellites [Eq. 6]
    ///   `phi(t) = phi(t_s) + SIGN(R, beta) * (t - t_s)`, `R = 0.20 deg/s`.
    /// GPS Block IIR satellites are not affected by Earth's shadow, but
    /// undergo a short (~15 min) midnight-turn maneuver near the orbit
    /// midnight point, limited by the maximum hardware yaw rate `R`.
    /// @param t Evaluation epoch
    /// @param ts Midnight-turn start epoch
    /// @param phi_ts Nominal yaw angle at `ts` [rad]
    /// @param beta Sun elevation angle above the orbital plane [rad]
    static Real GpsIIRMidnightTurnYawAngle(Real t, Real ts, Real phi_ts, Real beta);

    /// @brief Noon-turn maneuver model for GPS Block IIR satellites     [Eq. 7]
    ///   `phi(t) = phi(t_s) - SIGN(R, beta) * (t - t_s)`, `R = 0.20 deg/s`.
    /// @param t Evaluation epoch
    /// @param ts Noon-turn start epoch
    /// @param phi_ts Nominal yaw angle at `ts` [rad]
    /// @param beta Sun elevation angle above the orbital plane [rad]
    static Real GpsIIRNoonTurnYawAngle(Real t, Real ts, Real phi_ts, Real beta);

    // ---- GPS Block III -------------------------------------------------------
    // (Montenbruck et al., 2026)
    //
    // No manufacturer-published yaw-steering model exists for GPS Block III;
    // IGS analysis centers currently substitute the nominal, GPS IIR, or
    // GPS IIF laws, none of which match the attitude quaternions distributed
    // by JPL (believed to reflect the true on-board model) particularly well.
    // Montenbruck et al. (2026) instead adapt the Galileo-IOV-type "smoothed
    // yaw steering" law (cf. `GalileoIovEclipseYawAngle`) to GPS III's wider
    // collinearity zone and lower peak yaw rate (~0.09 deg/s, vs. ~0.20 deg/s
    // for IIR and ~0.11 deg/s for IIF), finding it to give the closest overall
    // approximation (RMS yaw-angle error 1.9-3.5 deg) of the JPL quaternions
    // among the candidate models considered -- and recommend it as a
    // transitional standard for OD/PPP processing of GPS III satellites.

    /// @brief Sun direction unit vector in the orbital reference frame, using
    /// the (along-track, orbit-normal, Earth-direction) axis convention and
    /// the midnight-referenced orbit angle `mu` directly (no detour through
    /// the noon-referenced `eta`), as introduced for the GPS III attitude
    /// model of Montenbruck et al. (2026)                            [Eq. 11]
    ///   `s = (cos(beta)*sin(mu), -sin(beta), cos(beta)*cos(mu))^T`
    /// Algebraically, `s == -GalileoSunVector(OrbitNoonAngle(mu), beta)`
    /// component-wise (the two papers' orbital-frame conventions differ by an
    /// overall sign), but both parameterizations yield the *same* nominal yaw
    /// angle:
    ///   `ATAN2(s_y, s_x) == ATAN2(-S_y, -S_x) == NominalYawAngle(beta, mu)`
    /// (dividing both `ATAN2` arguments by `cos(beta) > 0` collapses either
    /// expression to `ATAN2(-tan(beta), sin(mu))`, i.e. Eq. (1)/Eq. (12)).
    /// @param beta Sun elevation angle above the orbital plane [rad]
    /// @param mu Orbit angle from the midnight point [rad]
    static Vec3 Gps3SunVector(Real beta, Real mu);

    /// @brief Improved, rate-limited yaw-steering law for GPS Block III
    /// satellites during eclipse-season noon/midnight-turn maneuvers, per the
    /// GPS-III-adapted Galileo-IOV-type "smoothed yaw steering" model of
    /// Montenbruck et al. (2026)                          [Eq. 11-13, 17-18]
    ///   `phi(t) = ATAN2(s_y*, s_x)`                                 [Eq. 13]
    /// with the modified Sun-vector `y`-component
    ///   `s_y* = 0.5*(1+g)*SIGN(1,phi_nom)*gamma_y + 0.5*(1-g)*s_y` [Eq. 18]
    /// and IOV-type cosine weighting
    ///   `g = cos(pi*|s_x| / gamma_x)`                               [Eq. 17]
    /// `SIGN(1, phi_nom)` stands in for `sign(s_{y,0})` -- the sign of the
    /// Sun-vector `y`-component upon entry into the rectangular collinearity
    /// region `|s_x| < gamma_x`, `|s_y| < gamma_y` -- which is essentially
    /// constant over the short maneuver and shares the sign of `phi_nom`
    /// (Eq. 1 shows `phi_nom` and `beta` -- and hence `s_y = -sin(beta)` --
    /// always have opposite signs, so `SIGN(1,phi_nom) = sign(-beta) =
    /// sign(s_y)`). The collinearity-zone half-widths are GPS-III-specific
    /// (substantially wider than Galileo IOV's `gamma_x = sin(15deg)`,
    /// `gamma_y = sin(2deg)`, reflecting GPS III's lower peak yaw rate):
    ///   `gamma_x = sin(15 deg)`, `gamma_y = sin(5.8 deg)`
    /// @param beta Sun elevation angle above the orbital plane [rad]
    /// @param mu Orbit angle from the midnight point [rad]
    /// @param phi_nom Nominal yaw angle at the current epoch [rad]
    ///   (`NominalYawAngle(beta, mu)`, used only for its sign)
    static Real Gps3EclipseYawAngle(Real beta, Real mu, Real phi_nom);

    // ---- Galileo IOV ---------------------------------------------------------
    // (EGSC, 2017)

    /// @brief Sun reference vector in the orbital reference frame      [Eq. 9]
    ///   `S = (sin(eta)*cos(beta), sin(beta), cos(eta)*cos(beta))^T`
    /// where `eta` is the geocentric angle between the satellite and the
    /// orbit noon point (see `OrbitNoonAngle`): `eta = OrbitNoonAngle(mu)`.
    /// @param eta Orbit angle from the noon point [rad]
    /// @param beta Sun elevation angle above the orbital plane [rad]
    static Vec3 GalileoSunVector(Real eta, Real beta);

    /// @brief Nominal yaw angle for Galileo IOV satellites             [Eq. 8]
    ///   `phi(t) = ATAN2(-S_y/sqrt(1-S_z^2), -S_x/sqrt(1-S_z^2))`
    /// Algebraically identical to Eq. (1)/`NominalYawAngle(beta, mu)` once
    /// `eta = OrbitNoonAngle(mu)` is substituted via Eq. (9): with
    /// `sin(eta) = -sin(mu)`, `-S_x = -sin(eta)*cos(beta) = sin(mu)*cos(beta)`
    /// and `-S_y = -sin(beta)`, so `ATAN2(-S_y, -S_x) = ATAN2(-sin(beta),
    /// sin(mu)*cos(beta)) = ATAN2(-tan(beta), sin(mu))`.
    /// @param eta Orbit angle from the noon point [rad]
    /// @param beta Sun elevation angle above the orbital plane [rad]
    static Real GalileoIovNominalYawAngle(Real eta, Real beta);

    /// @brief Modeled (precise) yaw angle for Galileo IOV satellites during
    /// Earth's-shadow / noon-turn maneuvers                       [Eq. 10-11]
    ///   `phi(t) = ATAN2(-S_y'/sqrt(1-S_z^2), -S_x/sqrt(1-S_z^2))`
    /// with the modified Sun-vector `y`-component
    ///   `S_y' = [sin(2deg)*SIGN(1,phi_nom) + S_y]/2`
    ///         `+ [sin(2deg)*SIGN(1,phi_nom) - S_y]/2 * cos(pi*|S_x|/sin(15deg))`
    /// Triggered when `|beta| < 4.1 deg` and the orbit angle `mu` is within
    /// ~10 deg of the midnight/noon point.
    /// @param eta Orbit angle from the noon point [rad]
    /// @param beta Sun elevation angle above the orbital plane [rad]
    /// @param phi_nom Nominal yaw angle at the current epoch [rad]
    ///   (`NominalYawAngle(beta, mu)`, used only for its sign)
    static Real GalileoIovEclipseYawAngle(Real eta, Real beta, Real phi_nom);

    // ---- Galileo FOC ---------------------------------------------------------
    // (EGSC, 2017)

    /// @brief Yaw-steering law for Galileo FOC satellites during yaw
    /// maneuvers (Earth's shadow / noon-turn)                        [Eq. 12]
    ///   `phi(t) = 90deg*SIGN(1,phi_s)`
    ///           `+ (phi_s - 90deg*SIGN(1,phi_s)) * cos(2*pi/5656 * (t - t_s))`
    /// `5656` [s] is the maximum yaw-maneuver time constant. Triggered when
    /// `|beta| < 4.1 deg` and the orbit angle `mu` is within ~10 deg of the
    /// midnight/noon point.
    /// @param t Evaluation epoch
    /// @param ts Maneuver start (eclipsing) epoch
    /// @param phi_s Nominal yaw angle at `ts` [rad]
    static Real GalileoFocYawAngle(Real t, Real ts, Real phi_s);

    // ---- BDS-3 CAST (WHU model) ----------------------------------------------
    // (Wang et al., 2018)

    /// @brief WHU yaw-steering model for BDS-3 CAST IGSO satellites    [Eq. 13]
    ///   Same functional form as `GalileoFocYawAngle` (Eq. 12) with maximum
    ///   yaw-maneuver time `5740` [s] (maximum hardware yaw rate
    ///   `~0.20 deg/s`, midnight/noon turn lasting up to ~96 min).
    /// @param t Evaluation epoch
    /// @param ts Maneuver start (eclipsing) epoch
    /// @param phi_s Nominal yaw angle at `ts` [rad]
    static Real Bds3CastIgsoYawAngle(Real t, Real ts, Real phi_s);

    /// @brief WHU yaw-steering model for BDS-3 CAST MEO satellites     [Eq. 14]
    ///   Same functional form as `GalileoFocYawAngle` (Eq. 12) with maximum
    ///   yaw-maneuver time `3090` [s] (midnight/noon turn lasting ~52 min).
    /// @param t Evaluation epoch
    /// @param ts Maneuver start (eclipsing) epoch
    /// @param phi_s Nominal yaw angle at `ts` [rad]
    static Real Bds3CastMeoYawAngle(Real t, Real ts, Real phi_s);

    // ---- BDS-3 SECM (CSNO / MCSNO models) ------------------------------------
    // (CSNO, 2019; Xia et al., 2018; Yang et al., 2023)

    /// @brief CSNO model for BDS-3 SECM MEO satellites near zero-beta    [Eq. 15]
    ///   `phi(t) = ATAN2(-tan(3deg), sin(mu))`  for `0 < beta <= 3 deg`
    ///   `phi(t) = ATAN2( tan(3deg), sin(mu))`  for `-3 deg <= beta < 0`
    /// (`beta == 0` is treated as the `beta < 0` branch; this discontinuity
    /// at `beta = 0` is the "reverse yaw maneuver near zero beta-angle"
    /// reported by Xie et al. (2022) and is what `Bds3SecmMcsnoYawAngle`
    /// (the modified CSNO / "MCSNO" model) corrects).
    /// @param beta Sun elevation angle above the orbital plane [rad]
    ///   (model is only meaningful for `|beta| <= 3 deg`)
    /// @param mu Orbit angle from the midnight point [rad]
    static Real Bds3SecmCsnoYawAngle(Real beta, Real mu);

    /// @brief Modified-CSNO ("MCSNO") model for BDS-3 SECM MEO satellites,
    /// which replaces the CSNO model's discontinuous sign flip at `beta = 0`
    /// with a smooth linear transition of the yaw angle           [Eq. 16]
    ///
    ///   `phi(t) = ATAN2(-tan(3deg), sin(mu))`              for `t < t0, beta > 0`
    ///   `phi(t) = ATAN2( tan(3deg), sin(mu))`              for `t < t0, beta < 0`
    ///   `phi(t) = ATAN2(tan(SIGN(3deg, beta_dot)), sin(mu))`     for `t0 <= t < ts`
    ///   `phi(t) = phi(ts) - SIGN(0.055deg/s, beta_dot)`
    ///             `* (mu - mu_ts) / mu_dot`                      for `ts <= t < te`
    ///   `phi(t) = ATAN2(-tan(SIGN(3deg, beta_dot)), sin(mu))`    for `t > te`
    ///
    /// where `t0` is the epoch at which `beta` crosses zero, `beta_dot` is
    /// the rate of `beta`, `mu_dot` is the rate of `mu`, and `[ts, te)` is
    /// the ~3-minute linear-transition window during which the yaw angle
    /// ramps at the fixed rate `0.055 deg/s` (Yang et al., 2023) -- the
    /// transition begins at `t0` if the orbit angle `mu` at `t0` exceeds
    /// `36.8 deg`, and at `ts` (with `mu(ts) ~ 36.8 deg`) otherwise.
    /// @param t Evaluation epoch
    /// @param t0 Epoch at which `beta = 0`
    /// @param ts Linear-transition start epoch
    /// @param te Linear-transition end epoch
    /// @param beta Sun elevation angle above the orbital plane [rad]
    /// @param beta_dot Rate of `beta` [rad/s]
    /// @param mu Orbit angle from the midnight point at `t` [rad]
    /// @param mu_ts Orbit angle from the midnight point at `ts` [rad]
    /// @param mu_dot Rate of `mu` (orbital angular rate) [rad/s]
    /// @param phi_ts Yaw angle at the transition-start epoch `ts` [rad]
    static Real Bds3SecmMcsnoYawAngle(Real t, Real t0, Real ts, Real te, Real beta, Real beta_dot,
                                      Real mu, Real mu_ts, Real mu_dot, Real phi_ts);

    // ---- FORTRAN-style helper -------------------------------------------------

    /// @brief FORTRAN `SIGN(A, B)` intrinsic: the magnitude of `a` with the
    /// sign of `b` (`+|a|` if `b >= 0`, else `-|a|`). Used throughout the
    /// reference's yaw-steering laws to select the maneuver direction from
    /// the sign of `beta` (or `phi_nom`, `phi_s`, `beta_dot`).
    /// @param a Value whose magnitude is used
    /// @param b Value whose sign is used
    /// @return  `+|a|` if `b >= 0`, else `-|a|`
    static Real Sign(Real a, Real b);

  private:
    /// @brief Shared functional form of Eqs. (12)-(14): a half-cosine
    /// transition of the yaw angle from `phi_s` to `90 deg * SIGN(1, phi_s)`
    /// and back, with period `period_s` (the maximum yaw-maneuver duration).
    ///
    /// Called by `GalileoFocYawAngle`, `Bds3CastIgsoYawAngle`, and
    /// `Bds3CastMeoYawAngle`, which differ only in `period_s`.
    /// @param t Evaluation epoch
    /// @param ts Maneuver start (eclipsing) epoch
    /// @param phi_s Nominal yaw angle at `ts` [rad]
    /// @param period_s Maximum yaw-maneuver duration [s]
    /// @return Yaw angle `phi(t)` [rad]
    static Real CosineTransitionYawAngle(Real t, Real ts, Real phi_s, double period_s);
  };

}  // namespace lupnt
