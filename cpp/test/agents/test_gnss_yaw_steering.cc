#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-9;

TEST_CASE("agents.gnss_yaw_steering.geometry") {
  // Representative GPS-like circular orbit & Sun position, in ECI [m]/[m/s].
  // Orbit in the x-y plane (orbit normal +z), satellite at the +x point of
  // the orbit, moving in +y -- a small inclination is added so beta != 0.
  Real r = 26560e3;
  Real v = std::sqrt(static_cast<double>(GM_EARTH) / r.val());
  Vec3 r_sat(r, 0.0, 0.0);
  Vec3 v_sat(0.0, v, 0.0);

  SECTION("BetaAngle is the elevation of the Sun above the orbital plane") {
    // Sun in the orbital plane (x-y): beta = 0
    Vec3 r_sun_in_plane(0.3 * AU, 0.4 * AU, 0.0);
    RequireNear(GnssYawSteering::BetaAngle(r_sat, v_sat, r_sun_in_plane), Real(0.0), epsilon);

    // Sun straight above the orbital plane (+z, same side as orbit normal r x v):
    // beta = +90 deg. Note: asin(x) is singular (infinite derivative) at
    // x = +/-1, so a ~1e-16 round-off in the dot product gets amplified to
    // ~sqrt(eps) ~ 1e-8 in the resulting angle -- a much looser tolerance is
    // needed for this exactly-at-the-pole geometry than for regular angles.
    const double pole_epsilon = 1e-7;
    Vec3 r_sun_above(0.0, 0.0, AU);
    RequireNear(GnssYawSteering::BetaAngle(r_sat, v_sat, r_sun_above), Real(PI_OVER_TWO),
                pole_epsilon);

    // Sun straight below the orbital plane: beta = -90 deg
    Vec3 r_sun_below(0.0, 0.0, -AU);
    RequireNear(GnssYawSteering::BetaAngle(r_sat, v_sat, r_sun_below), Real(-PI_OVER_TWO),
                pole_epsilon);
  }

  SECTION("OrbitAngle is zero at the midnight point and +/-pi at the noon point") {
    // Sun along -x: the midnight point (antipodal to the sub-solar point,
    // projected onto the orbital plane) is at +x -- exactly where r_sat is.
    Vec3 r_sun_midnight(-AU, 0.0, 0.0);
    RequireNear(GnssYawSteering::OrbitAngle(r_sat, v_sat, r_sun_midnight), Real(0.0), epsilon);

    // Sun along +x: the midnight point is at -x; the satellite (at +x) is at
    // the noon point, i.e. mu = +/- pi.
    Vec3 r_sun_noon(AU, 0.0, 0.0);
    Real mu_noon = GnssYawSteering::OrbitAngle(r_sat, v_sat, r_sun_noon);
    RequireNear(Real(abs(mu_noon)), Real(PI), epsilon);

    // Sun along +y: the midnight point is at -y, i.e. 90 deg (in the
    // direction of motion, +y) before the satellite's position (+x) -> mu = +90 deg
    Vec3 r_sun_dawn(0.0, AU, 0.0);
    RequireNear(GnssYawSteering::OrbitAngle(r_sat, v_sat, r_sun_dawn), Real(PI_OVER_TWO), epsilon);
  }

  SECTION("OrbitNoonAngle converts a midnight- to a noon-referenced orbit angle") {
    RequireNear(GnssYawSteering::OrbitNoonAngle(Real(0.0)), Real(-PI), epsilon);
    RequireNear(GnssYawSteering::OrbitNoonAngle(Real(PI_OVER_TWO)), Real(-PI_OVER_TWO), epsilon);
    RequireNear(GnssYawSteering::OrbitNoonAngle(Real(PI)), Real(0.0), epsilon);
  }
}

TEST_CASE("agents.gnss_yaw_steering.nominal_model") {
  // Eq. (1): the sign of the nominal yaw angle is always opposite that of beta
  Real mu = 30.0 * RAD;
  Real beta_pos = 5.0 * RAD;
  Real beta_neg = -5.0 * RAD;

  Real phi_pos_beta = GnssYawSteering::NominalYawAngle(beta_pos, mu);
  Real phi_neg_beta = GnssYawSteering::NominalYawAngle(beta_neg, mu);
  REQUIRE(phi_pos_beta.val() < 0.0);
  REQUIRE(phi_neg_beta.val() > 0.0);
  // The nominal model is symmetric in beta -> -beta (Eq. 1: -tan(-b) = tan(b))
  RequireNear(phi_pos_beta, -phi_neg_beta, epsilon);

  // At the orbit midnight (mu = 0) and noon (mu = +/- pi) points, sin(mu) = 0
  // and the nominal yaw is +/- 90 deg (a vertical jump in the nominal model,
  // which is what triggers the midnight/noon-turn maneuvers in Section 2).
  RequireNear(Real(abs(GnssYawSteering::NominalYawAngle(beta_pos, Real(0.0)))), Real(PI_OVER_TWO),
              epsilon);
  RequireNear(Real(abs(GnssYawSteering::NominalYawAngle(beta_pos, Real(PI)))), Real(PI_OVER_TWO),
              epsilon);

  // Eq. (2): the nominal yaw rate is finite away from beta = 0 / mu = 0, pi
  Real mu_dot = TWO_PI / 43200.0;  // ~ MEO orbital rate [rad/s]
  Real phi_dot = GnssYawSteering::NominalYawRate(beta_pos, mu, mu_dot);
  REQUIRE(std::isfinite(phi_dot.val()));

  // Numerically differentiate the nominal yaw angle and compare to Eq. (2)
  double dmu = 1e-6;
  Real phi_p = GnssYawSteering::NominalYawAngle(beta_pos, mu + dmu);
  Real phi_m = GnssYawSteering::NominalYawAngle(beta_pos, mu - dmu);
  Real dphi_dmu = (phi_p - phi_m) / (2.0 * dmu);
  RequireNear(phi_dot, dphi_dmu * mu_dot, 1e-6);
}

TEST_CASE("agents.gnss_yaw_steering.gps") {
  Real ts = 1000.0;
  Real te = 4000.0;
  Real phi_ts = 170.0 * RAD;
  Real phi_te = -170.0 * RAD;

  SECTION("GPS Block IIF shadow-crossing model interpolates linearly  [Eq. 3-4]") {
    RequireNear(GnssYawSteering::GpsIIFShadowYawAngle(ts, ts, te, phi_ts, phi_te), phi_ts, epsilon);
    RequireNear(GnssYawSteering::GpsIIFShadowYawAngle(te, ts, te, phi_ts, phi_te), phi_te, epsilon);
    Real t_mid = (ts + te) / 2.0;
    RequireNear(GnssYawSteering::GpsIIFShadowYawAngle(t_mid, ts, te, phi_ts, phi_te),
                (phi_ts + phi_te) / 2.0, epsilon);
  }

  SECTION("GPS Block IIF noon-turn model uses the maximum hardware yaw rate  [Eq. 5]") {
    Real beta = 1.0 * RAD;
    Real t = ts + 100.0;
    Real phi = GnssYawSteering::GpsIIFNoonTurnYawAngle(t, ts, phi_ts, beta);
    constexpr double kR = 0.11 * RAD;
    // beta + 0.7 deg > 0 -> SIGN(R, .) = +R -> phi decreases at rate R
    RequireNear(phi, phi_ts - Real(kR) * (t - ts), epsilon);
    REQUIRE(phi.val() < phi_ts.val());

    // Flipping the sign of (beta + 0.7 deg) flips the maneuver direction
    Real beta_neg = -2.0 * RAD;  // beta + 0.7 deg < 0
    Real phi_flipped = GnssYawSteering::GpsIIFNoonTurnYawAngle(t, ts, phi_ts, beta_neg);
    RequireNear(phi_flipped, phi_ts + Real(kR) * (t - ts), epsilon);
  }

  SECTION("GPS Block IIR midnight/noon-turn models use a slower hardware yaw rate  [Eq. 6-7]") {
    Real beta = 0.3 * RAD;
    Real t = ts + 50.0;
    constexpr double kR = 0.20 * RAD;

    Real phi_midnight = GnssYawSteering::GpsIIRMidnightTurnYawAngle(t, ts, phi_ts, beta);
    RequireNear(phi_midnight, phi_ts + Real(kR) * (t - ts), epsilon);

    Real phi_noon = GnssYawSteering::GpsIIRNoonTurnYawAngle(t, ts, phi_ts, beta);
    RequireNear(phi_noon, phi_ts - Real(kR) * (t - ts), epsilon);

    // beta < 0 flips the maneuver direction (SIGN(R, beta))
    Real beta_neg = -0.3 * RAD;
    RequireNear(GnssYawSteering::GpsIIRMidnightTurnYawAngle(t, ts, phi_ts, beta_neg),
                phi_ts - Real(kR) * (t - ts), epsilon);
  }
}

TEST_CASE("agents.gnss_yaw_steering.galileo") {
  SECTION("GalileoSunVector is a unit vector and matches Eq. (9) component-wise") {
    Real eta = 20.0 * RAD;
    Real beta = 3.0 * RAD;
    Vec3 S = GnssYawSteering::GalileoSunVector(eta, beta);
    RequireNear(Real(S.norm()), Real(1.0), 1e-12);
    RequireNear(S.x(), sin(eta) * cos(beta), epsilon);
    RequireNear(S.y(), sin(beta), epsilon);
    RequireNear(S.z(), cos(eta) * cos(beta), epsilon);
  }

  SECTION("Galileo IOV nominal yaw model is consistent with the generic nominal model") {
    // Eq. (8)/(9) are *algebraically identical* to the generic nominal model
    // phi = atan2(-tan(beta), sin(mu))  (Eq. 1) once eta = OrbitNoonAngle(mu)
    // = mu - pi is substituted:
    //   -S_y = -sin(beta)
    //   -S_x = -sin(eta)*cos(beta) = -sin(mu - pi)*cos(beta) = sin(mu)*cos(beta)
    // so (dropping the common positive scale factor 1/sqrt(1 - S_z^2)):
    //   ATAN2(-S_y, -S_x) = ATAN2(-sin(beta), sin(mu)*cos(beta))
    //                     = ATAN2(-sin(beta)/cos(beta), sin(mu))   [cos(beta) > 0]
    //                     = ATAN2(-tan(beta), sin(mu))  =  Eq. (1)
    // The two expressions therefore match exactly (not merely up to an
    // additive offset), which is the "equivalent to Eq. (1)" property the
    // paper states right after Eq. (9).
    Real beta = 1.0 * RAD;
    Real mu = 50.0 * RAD;
    Real eta = GnssYawSteering::OrbitNoonAngle(mu);

    Real phi_galileo = GnssYawSteering::GalileoIovNominalYawAngle(eta, beta);
    Real phi_nominal = GnssYawSteering::NominalYawAngle(beta, mu);
    RequireNear(phi_galileo, phi_nominal, 1e-9);
  }

  SECTION("Galileo IOV eclipse model produces a finite, valid yaw angle  [Eq. 10-11]") {
    Real beta = 0.5 * RAD;
    Real mu = 5.0 * RAD;
    Real eta = GnssYawSteering::OrbitNoonAngle(mu);
    Real phi_nom = GnssYawSteering::NominalYawAngle(beta, mu);

    Real phi = GnssYawSteering::GalileoIovEclipseYawAngle(eta, beta, phi_nom);
    REQUIRE(std::isfinite(phi.val()));
    REQUIRE(phi.val() >= -PI);
    REQUIRE(phi.val() <= PI);
  }

  SECTION(
      "Galileo FOC yaw-steering law oscillates between phi_s and 90 deg * SIGN(1, phi_s)  "
      "[Eq. 12]") {
    Real ts = 0.0;
    Real phi_s = 150.0 * RAD;
    constexpr double kPeriod = 5656.0;

    // At t = ts: cos(0) = 1 -> phi = phi_s
    RequireNear(GnssYawSteering::GalileoFocYawAngle(ts, ts, phi_s), phi_s, epsilon);

    // At t = ts + period/2: cos(pi) = -1 -> phi = 2*90deg*SIGN(1,phi_s) - phi_s
    Real t_half = ts + kPeriod / 2.0;
    Real expected_half = Real(PI) * GnssYawSteering::Sign(1.0, phi_s) - phi_s;
    RequireNear(GnssYawSteering::GalileoFocYawAngle(t_half, ts, phi_s), expected_half, epsilon);

    // At t = ts + period: back to phi_s (full cosine cycle)
    RequireNear(GnssYawSteering::GalileoFocYawAngle(ts + kPeriod, ts, phi_s), phi_s, epsilon);

    // phi_s > 0 -> SIGN(1, phi_s) = +1 -> the law oscillates toward +90 deg
    REQUIRE(GnssYawSteering::Sign(1.0, phi_s).val() == 1.0);
  }
}

TEST_CASE("agents.gnss_yaw_steering.gps3") {
  SECTION("Gps3SunVector is a unit vector and matches Eq. (11) component-wise") {
    Real beta = 3.0 * RAD;
    Real mu = 20.0 * RAD;
    Vec3 s = GnssYawSteering::Gps3SunVector(beta, mu);
    RequireNear(Real(s.norm()), Real(1.0), 1e-12);
    RequireNear(s.x(), cos(beta) * sin(mu), epsilon);
    RequireNear(s.y(), -sin(beta), epsilon);
    RequireNear(s.z(), cos(beta) * cos(mu), epsilon);
  }

  SECTION("Gps3SunVector nominal yaw (Eq. 12) is consistent with the generic nominal model") {
    // Eq. (12): phi_nom = ATAN2(s_y, s_x). Substituting Eq. (11):
    //   ATAN2(s_y, s_x) = ATAN2(-sin(beta), cos(beta)*sin(mu))
    //                   = ATAN2(-sin(beta)/cos(beta), sin(mu))   [cos(beta) > 0]
    //                   = ATAN2(-tan(beta), sin(mu))  =  Eq. (1) = NominalYawAngle
    // i.e. the new paper's direct-mu Sun-vector parameterization yields
    // *exactly* the same nominal yaw angle as the generic nominal model (and
    // as the eta-based Galileo parameterization, since s == -GalileoSunVector
    // component-wise but the sign cancels in the ATAN2 ratio).
    Real beta = 1.0 * RAD;
    Real mu = 50.0 * RAD;
    Vec3 s = GnssYawSteering::Gps3SunVector(beta, mu);

    Real phi_gps3 = atan2(s.y(), s.x());
    Real phi_nominal = GnssYawSteering::NominalYawAngle(beta, mu);
    RequireNear(phi_gps3, phi_nominal, 1e-9);
  }

  SECTION("Gps3EclipseYawAngle produces a finite, valid yaw angle  [Eq. 13, 17-18]") {
    Real beta = 0.5 * RAD;
    Real mu = 5.0 * RAD;
    Real phi_nom = GnssYawSteering::NominalYawAngle(beta, mu);

    Real phi = GnssYawSteering::Gps3EclipseYawAngle(beta, mu, phi_nom);
    REQUIRE(std::isfinite(phi.val()));
    REQUIRE(phi.val() >= -PI);
    REQUIRE(phi.val() <= PI);
  }

  SECTION(
      "Gps3EclipseYawAngle reduces to known limiting forms at the collinearity-region "
      "boundaries  [Eq. 17-18]") {
    // At the collinearity-zone edge |s_x| = gamma_x = sin(15 deg): g = cos(pi) = -1,
    // so Eq. (18) collapses to s_y* = s_y, i.e. the eclipse model reduces to the
    // (unmodified) nominal Sun-vector and hence to the nominal yaw angle.
    Real gamma_x = sind(15.0);
    Real beta = 2.0 * RAD;
    // Choose mu so that cos(beta)*sin(mu) = +gamma_x  =>  sin(mu) = gamma_x / cos(beta)
    Real mu = asin(gamma_x / cos(beta)).val();
    Real phi_nom = GnssYawSteering::NominalYawAngle(beta, mu);

    Vec3 s = GnssYawSteering::Gps3SunVector(beta, mu);
    RequireNear(Real(abs(s.x())), gamma_x, 1e-9);

    Real phi_edge = GnssYawSteering::Gps3EclipseYawAngle(beta, mu, phi_nom);
    Real phi_expected = atan2(s.y(), s.x());
    RequireNear(phi_edge, phi_expected, 1e-9);
    RequireNear(phi_edge, phi_nom, 1e-9);

    // At the orbit noon/midnight point s_x = 0: g = cos(0) = 1, so Eq. (18)
    // collapses to s_y* = SIGN(1, phi_nom) * gamma_y, the full collinearity-
    // region clamp (the maximal rate-limiting departure from the nominal law).
    Real gamma_y = sind(5.8);
    Real mu0 = 0.0;  // s_x = cos(beta)*sin(0) = 0
    Real phi_nom0 = GnssYawSteering::NominalYawAngle(beta, mu0);
    Real phi0 = GnssYawSteering::Gps3EclipseYawAngle(beta, mu0, phi_nom0);
    Real sgn0 = GnssYawSteering::Sign(1.0, phi_nom0);
    Real phi0_expected = atan2(sgn0 * gamma_y, Real(0.0));
    RequireNear(phi0, phi0_expected, 1e-9);
  }
}

TEST_CASE("agents.gnss_yaw_steering.bds3_cast_whu") {
  // Eqs. (13)-(14) share the functional form of Eq. (12), only the maximum
  // yaw-maneuver-time constant differs (5740 s for IGSO, 3090 s for MEO).
  Real ts = 200.0;
  Real phi_s = -160.0 * RAD;

  Real phi_igso_start = GnssYawSteering::Bds3CastIgsoYawAngle(ts, ts, phi_s);
  Real phi_meo_start = GnssYawSteering::Bds3CastMeoYawAngle(ts, ts, phi_s);
  RequireNear(phi_igso_start, phi_s, epsilon);
  RequireNear(phi_meo_start, phi_s, epsilon);

  // The two models agree at t = ts (cos(0) = 1) but differ at other epochs
  // because of the different oscillation periods.
  Real t = ts + 1000.0;
  Real phi_igso = GnssYawSteering::Bds3CastIgsoYawAngle(t, ts, phi_s);
  Real phi_meo = GnssYawSteering::Bds3CastMeoYawAngle(t, ts, phi_s);
  REQUIRE(std::abs(phi_igso.val() - phi_meo.val()) > 1e-3);

  // Both models complete a full oscillation cycle and return to phi_s after
  // their respective periods.
  RequireNear(GnssYawSteering::Bds3CastIgsoYawAngle(ts + 5740.0, ts, phi_s), phi_s, epsilon);
  RequireNear(GnssYawSteering::Bds3CastMeoYawAngle(ts + 3090.0, ts, phi_s), phi_s, epsilon);
}

TEST_CASE("agents.gnss_yaw_steering.bds3_secm") {
  SECTION("CSNO model: the yaw target sign is opposite that of beta near zero-beta  [Eq. 15]") {
    Real mu = 10.0 * RAD;
    constexpr double kThreeDeg = 3.0 * RAD;

    Real beta_pos = 1.0 * RAD;
    Real phi_pos = GnssYawSteering::Bds3SecmCsnoYawAngle(beta_pos, mu);
    RequireNear(phi_pos, atan2(Real(-tan(kThreeDeg)), sin(mu)), epsilon);

    Real beta_neg = -1.0 * RAD;
    Real phi_neg = GnssYawSteering::Bds3SecmCsnoYawAngle(beta_neg, mu);
    RequireNear(phi_neg, atan2(Real(tan(kThreeDeg)), sin(mu)), epsilon);

    // The discontinuous sign flip across beta = 0 is exactly the "reverse
    // direction of yaw maneuver" reported by Xie et al. (2022) -- the two
    // branches give different (non-equal) yaw targets.
    REQUIRE(std::abs(phi_pos.val() - phi_neg.val()) > 1e-3);
  }

  SECTION(
      "MCSNO model matches CSNO before t0 and is continuous through the linear "
      "transition  [Eq. 16]") {
    Real t0 = 1000.0;
    Real ts = 1200.0;
    Real te = 1380.0;  // ~3 minute transition window
    Real mu_dot = TWO_PI / 43200.0;
    constexpr double kThreeDeg = 3.0 * RAD;
    constexpr double kTransitionRate = 0.055 * RAD;

    // --- Before t0: matches the (discontinuous) CSNO branches directly
    Real beta_before_pos = 0.5 * RAD;
    Real mu_before = 100.0 * RAD;
    Real phi_before_pos = GnssYawSteering::Bds3SecmMcsnoYawAngle(
        Real(500.0), t0, ts, te, beta_before_pos, Real(1e-7), mu_before, Real(0.0), mu_dot,
        Real(0.0));
    RequireNear(phi_before_pos, atan2(Real(-tan(kThreeDeg)), sin(mu_before)), epsilon);

    Real beta_before_neg = -0.5 * RAD;
    Real phi_before_neg = GnssYawSteering::Bds3SecmMcsnoYawAngle(
        Real(500.0), t0, ts, te, beta_before_neg, Real(-1e-7), mu_before, Real(0.0), mu_dot,
        Real(0.0));
    RequireNear(phi_before_neg, atan2(Real(tan(kThreeDeg)), sin(mu_before)), epsilon);

    // --- [t0, ts): target yaw angle set by SIGN(3 deg, beta_dot)
    Real beta_dot = -2e-7;  // beta decreasing through zero -> SIGN(3 deg, beta_dot) = -3 deg
    Real mu_t0_ts = 80.0 * RAD;
    Real phi_t0_ts
        = GnssYawSteering::Bds3SecmMcsnoYawAngle(Real(1100.0), t0, ts, te, Real(0.05 * RAD),
                                                 beta_dot, mu_t0_ts, Real(0.0), mu_dot, Real(0.0));
    RequireNear(phi_t0_ts,
                atan2(Real(tan(GnssYawSteering::Sign(kThreeDeg, beta_dot))), sin(mu_t0_ts)),
                epsilon);

    // --- [ts, te): linear transition continuous with the [t0, ts) branch at t = ts
    //
    // Per Eq. (16), the function's `phi_ts` parameter is "the nominal yaw
    // angle at the [transition] starting time ts" -- i.e. it is *supplied* by
    // the caller (typically the [t0, ts) branch's value evaluated at mu =
    // mu_ts), not computed internally. Note that calling the function at
    // exactly `t = ts` lands in the `[ts, te)` branch (the `t < ts` test is
    // false at equality) and simply returns the supplied `phi_ts` unchanged,
    // since the linear-transition term `(mu - mu_ts)` vanishes there -- so we
    // compute the continuity-enforcing boundary value directly and supply it,
    // mirroring how an upstream caller would chain the two branches together.
    Real mu_ts = 90.0 * RAD;
    Real expected_phi_ts = atan2(Real(tan(GnssYawSteering::Sign(kThreeDeg, beta_dot))), sin(mu_ts));

    Real phi_ts = GnssYawSteering::Bds3SecmMcsnoYawAngle(ts, t0, ts, te, Real(0.0 * RAD), beta_dot,
                                                         mu_ts, mu_ts, mu_dot, expected_phi_ts);
    RequireNear(phi_ts, expected_phi_ts, epsilon);

    // Linear transition rate matches the documented 0.055 deg/s
    Real mu_t = mu_ts + mu_dot * 60.0;  // ~ orbit angle 60 s after ts
    Real phi_t = GnssYawSteering::Bds3SecmMcsnoYawAngle(
        ts + 60.0, t0, ts, te, Real(0.0 * RAD), beta_dot, mu_t, mu_ts, mu_dot, expected_phi_ts);
    RequireNear(phi_t, expected_phi_ts - GnssYawSteering::Sign(kTransitionRate, beta_dot) * 60.0,
                1e-9);

    // --- After te: target yaw angle flips sign relative to the [t0, ts) branch
    Real mu_after = 95.0 * RAD;
    Real phi_after = GnssYawSteering::Bds3SecmMcsnoYawAngle(
        te + 10.0, t0, ts, te, Real(0.0 * RAD), beta_dot, mu_after, mu_ts, mu_dot, Real(0.0));
    RequireNear(phi_after,
                atan2(Real(-tan(GnssYawSteering::Sign(kThreeDeg, beta_dot))), sin(mu_after)),
                epsilon);
  }
}

TEST_CASE("agents.gnss_yaw_steering.sign") {
  // FORTRAN SIGN(a, b): magnitude of `a` with the sign of `b`
  RequireNear(GnssYawSteering::Sign(Real(2.0), Real(3.0)), Real(2.0), epsilon);
  RequireNear(GnssYawSteering::Sign(Real(2.0), Real(-3.0)), Real(-2.0), epsilon);
  RequireNear(GnssYawSteering::Sign(Real(-2.0), Real(3.0)), Real(2.0), epsilon);
  RequireNear(GnssYawSteering::Sign(Real(-2.0), Real(-3.0)), Real(-2.0), epsilon);
  // b = 0 is treated as non-negative (FORTRAN convention: SIGN(a, 0) = |a|)
  RequireNear(GnssYawSteering::Sign(Real(-2.0), Real(0.0)), Real(2.0), epsilon);
}
