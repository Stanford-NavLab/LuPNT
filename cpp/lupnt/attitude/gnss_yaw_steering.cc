/**
 * @file gnss_yaw_steering.cc
 * @author Stanford NAV LAB
 * @brief GNSS satellite yaw-attitude steering laws
 * @version 0.1
 * @date 2026-06-07
 *
 * @copyright Copyright (c) 2026
 *
 * See `gnss_yaw_steering.h` for the reference (Cheng et al., 2025) and the
 * equation numbers each function mirrors.
 */
#include "lupnt/attitude/gnss_yaw_steering.h"

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {

  namespace {
    // Empirical constants from the reference, converted from degrees (and
    // deg/s) to the radian-based convention used throughout `lupnt`.
    constexpr double kNinetyDeg = PI_OVER_TWO;            // 90 deg
    constexpr double kThreeDeg = 3.0 * RAD;               // 3 deg          (Eq. 15-16)
    constexpr double kGpsIIFNoonTurnRate = 0.11 * RAD;    // 0.11 deg/s     (Eq. 5)
    constexpr double kGpsIIFBetaBias = 0.7 * RAD;         // 0.7 deg        (Eq. 5)
    constexpr double kGpsIIRYawRate = 0.20 * RAD;         // 0.20 deg/s     (Eq. 6-7)
    constexpr double kMcsnoTransitionRate = 0.055 * RAD;  // 0.055 deg/s    (Eq. 16)
    constexpr double kGalileoFocPeriod = 5656.0;          // [s]            (Eq. 12)
    constexpr double kBds3IgsoPeriod = 5740.0;            // [s]            (Eq. 13)
    constexpr double kBds3MeoPeriod = 3090.0;             // [s]            (Eq. 14)
  }  // namespace

  // ---- Satellite-Sun-Earth geometry ----------------------------------------

  Real GnssYawSteering::BetaAngle(const Vec3& r_sat, const Vec3& v_sat, const Vec3& r_sun) {
    Vec3 n = r_sat.cross(v_sat).normalized();  // Orbit-normal (angular-momentum) direction
    Vec3 s_hat = r_sun.normalized();
    return safe_asin(n.dot(s_hat));
  }

  Real GnssYawSteering::OrbitAngle(const Vec3& r_sat, const Vec3& v_sat, const Vec3& r_sun) {
    Vec3 n = r_sat.cross(v_sat).normalized();
    Vec3 s_hat = r_sun.normalized();
    Vec3 r_hat = r_sat.normalized();
    // Anti-Sun direction projected onto the orbital plane: points toward the
    // orbit "midnight point" (the point of the orbit directly behind Earth,
    // as seen from the Sun).
    Vec3 m = -s_hat - (-s_hat).dot(n) * n;
    Vec3 m_hat = m.normalized();
    // Signed angle from `m_hat` to `r_hat`, about the orbit-normal axis `n`
    // (i.e. positive in the direction of orbital motion).
    return atan2(m_hat.cross(r_hat).dot(n), m_hat.dot(r_hat));
  }

  Real GnssYawSteering::OrbitNoonAngle(Real mu) { return WrapToPi(mu - PI); }

  // ---- Nominal yaw-attitude model (Kouba, 2009), Eq. (1)-(2) ---------------

  Real GnssYawSteering::NominalYawAngle(Real beta, Real mu) { return atan2(-tan(beta), sin(mu)); }

  Real GnssYawSteering::NominalYawRate(Real beta, Real mu, Real mu_dot) {
    Real tan_beta = tan(beta);
    Real sin_mu = sin(mu);
    return mu_dot * tan_beta * cos(mu) / (sin_mu * sin_mu + tan_beta * tan_beta);
  }

  // ---- GPS Block IIF --------------------------------------------------------

  Real GnssYawSteering::GpsIIFShadowYawAngle(Real t, Real ts, Real te, Real phi_ts, Real phi_te) {
    Real phi_dot_star = (phi_te - phi_ts) / (te - ts);  // Eq. (3)
    return phi_ts + phi_dot_star * (t - ts);            // Eq. (4)
  }

  Real GnssYawSteering::GpsIIFNoonTurnYawAngle(Real t, Real ts, Real phi_ts, Real beta) {
    return phi_ts - Sign(kGpsIIFNoonTurnRate, beta + kGpsIIFBetaBias) * (t - ts);  // Eq. (5)
  }

  // ---- GPS Block IIR --------------------------------------------------------

  Real GnssYawSteering::GpsIIRMidnightTurnYawAngle(Real t, Real ts, Real phi_ts, Real beta) {
    return phi_ts + Sign(kGpsIIRYawRate, beta) * (t - ts);  // Eq. (6)
  }

  Real GnssYawSteering::GpsIIRNoonTurnYawAngle(Real t, Real ts, Real phi_ts, Real beta) {
    return phi_ts - Sign(kGpsIIRYawRate, beta) * (t - ts);  // Eq. (7)
  }

  // ---- GPS Block III (Montenbruck et al., 2026) ----------------------------

  Vec3 GnssYawSteering::Gps3SunVector(Real beta, Real mu) {
    Real cos_beta = cos(beta);
    return Vec3(cos_beta * sin(mu), -sin(beta), cos_beta * cos(mu));  // Eq. (11)
  }

  Real GnssYawSteering::Gps3EclipseYawAngle(Real beta, Real mu, Real phi_nom) {
    Vec3 s = Gps3SunVector(beta, mu);
    Real sx = s.x(), sy = s.y();

    // GPS-III-adapted Galileo-IOV-type collinearity-region half-widths
    // (Montenbruck et al., 2026): +/- 15 deg in orbit angle and +/- 5.8 deg
    // in Sun elevation -- substantially wider than Galileo IOV's
    // (+/- 15 deg, +/- 2 deg), reflecting GPS III's lower peak yaw rate.
    Real sin_15deg = sind(15.0);
    Real sin_5p8deg = sind(5.8);

    // The rate-limited law is defined ONLY inside the collinearity region; outside it
    // (|s_x| >= gamma_x, i.e. away from the noon/midnight turn, or |beta| >= beta_0) the
    // satellite follows nominal yaw steering. Evaluated continuously, so this gate is
    // what keeps the yaw nominal away from the turns.
    if (abs(sx) >= sin_15deg || abs(beta) >= 5.8 * RAD) return phi_nom;

    Real sgn = Sign(1.0, phi_nom);
    Real g = cos(PI * abs(sx) / sin_15deg);  // Eq. (17)
    // Eq. (18): IOV-type modified Sun-vector y-component
    Real sy_p = 0.5 * (1.0 + g) * sgn * sin_5p8deg + 0.5 * (1.0 - g) * sy;
    return atan2(sy_p, sx);  // Eq. (13)
  }

  // ---- Galileo IOV ----------------------------------------------------------

  Vec3 GnssYawSteering::GalileoSunVector(Real eta, Real beta) {
    Real cos_beta = cos(beta);
    return Vec3(sin(eta) * cos_beta, sin(beta), cos(eta) * cos_beta);  // Eq. (9)
  }

  Real GnssYawSteering::GalileoIovNominalYawAngle(Real eta, Real beta) {
    Vec3 S = GalileoSunVector(eta, beta);
    Real denom = sqrt(1.0 - S.z() * S.z());
    return atan2(-S.y() / denom, -S.x() / denom);  // Eq. (8)
  }

  Real GnssYawSteering::GalileoIovEclipseYawAngle(Real eta, Real beta, Real phi_nom) {
    Vec3 S = GalileoSunVector(eta, beta);
    Real Sx = S.x(), Sy = S.y(), Sz = S.z();

    Real sin_2deg = sind(2.0);
    Real sin_15deg = sind(15.0);
    // Defined ONLY inside the collinearity region (Galileo IOV: |S_x| < gamma_x = sin 15 deg
    // and |beta| < gamma_y = 2 deg); outside it the satellite follows nominal yaw steering.
    if (abs(Sx) >= sin_15deg || abs(beta) >= 2.0 * RAD) return phi_nom;

    Real sgn = Sign(1.0, phi_nom);
    // Eq. (11): modified Sun-vector y-component. The target term carries a minus sign
    // relative to Cheng et al. Eq. (11) to match lupnt's Sun-vector convention, which is
    // negated w.r.t. the paper so the nominal law (`GalileoIovNominalYawAngle`) reduces
    // exactly to Kouba Eq. (1). With this sign the rate-limited slew follows the nominal
    // hemisphere (short-way turn), matching the paper's Fig. 2; algebraically this whole
    // branch is identical to the GPS-III/Montenbruck IOV form `atan2(s*_y, s_x)` with
    // gamma_y = sin 2 deg (verified numerically).
    Real Sy_p = (-sin_2deg * sgn + Sy) / 2.0
                + (-sin_2deg * sgn - Sy) / 2.0 * cos(PI * abs(Sx) / sin_15deg);

    Real denom = sqrt(1.0 - Sz * Sz);
    return atan2(-Sy_p / denom, -Sx / denom);  // Eq. (10)
  }

  // ---- Galileo FOC / BDS-3 CAST (WHU model): shared cosine-transition form -

  Real GnssYawSteering::CosineTransitionYawAngle(Real t, Real ts, Real phi_s, double period_s) {
    Real sgn = Sign(1.0, phi_s);
    Real phi_max = kNinetyDeg * sgn;
    return phi_max + (phi_s - phi_max) * cos(TWO_PI / period_s * (t - ts));
  }

  Real GnssYawSteering::GalileoFocYawAngle(Real t, Real ts, Real phi_s) {
    return CosineTransitionYawAngle(t, ts, phi_s, kGalileoFocPeriod);  // Eq. (12)
  }

  Real GnssYawSteering::Bds3CastIgsoYawAngle(Real t, Real ts, Real phi_s) {
    return CosineTransitionYawAngle(t, ts, phi_s, kBds3IgsoPeriod);  // Eq. (13)
  }

  Real GnssYawSteering::Bds3CastMeoYawAngle(Real t, Real ts, Real phi_s) {
    return CosineTransitionYawAngle(t, ts, phi_s, kBds3MeoPeriod);  // Eq. (14)
  }

  // ---- BDS-3 SECM (CSNO / MCSNO models) -------------------------------------

  Real GnssYawSteering::Bds3SecmCsnoYawAngle(Real beta, Real mu) {
    // Eq. (15): the threshold angle 3 deg sets the target yaw angle during
    // the yaw maneuver; the sign flips discontinuously across `beta = 0`
    // (the "reverse direction of yaw maneuver" reported by Xie et al., 2022,
    // and corrected by the MCSNO model, `Bds3SecmMcsnoYawAngle`, below).
    if (beta > 0.0) {
      return atan2(-tan(kThreeDeg), sin(mu));
    } else {
      return atan2(tan(kThreeDeg), sin(mu));
    }
  }

  Real GnssYawSteering::Bds3SecmMcsnoYawAngle(Real t, Real t0, Real ts, Real te, Real beta,
                                              Real beta_dot, Real mu, Real mu_ts, Real mu_dot,
                                              Real phi_ts) {
    // Eq. (16)
    if (t < t0) {
      return (beta > 0.0) ? atan2(-tan(kThreeDeg), sin(mu)) : atan2(tan(kThreeDeg), sin(mu));
    }

    Real tan_target = tan(Sign(kThreeDeg, beta_dot));
    if (t < ts) {
      return atan2(tan_target, sin(mu));
    } else if (t < te) {
      // Linear transition of the yaw angle at the fixed rate 0.055 deg/s;
      // `(mu - mu_ts) / mu_dot` converts the orbit-angle change since `ts`
      // back into an elapsed time, `~= (t - ts)`.
      return phi_ts - Sign(kMcsnoTransitionRate, beta_dot) * (mu - mu_ts) / mu_dot;
    } else {
      return atan2(-tan_target, sin(mu));
    }
  }

  // ---- FORTRAN-style helper --------------------------------------------------

  Real GnssYawSteering::Sign(Real a, Real b) { return (b >= 0.0) ? abs(a) : -abs(a); }

}  // namespace lupnt
