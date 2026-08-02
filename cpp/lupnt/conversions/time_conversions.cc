#include "lupnt/conversions/time_conversions.h"

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "lupnt/conversions/epoch.h"
#include "lupnt/core/error.h"
#include "lupnt/environment/body.h"
#include "lupnt/environment/solar_system.h"
#include "lupnt/interfaces/eop.h"
#include "lupnt/interfaces/kernels.h"
#include "lupnt/interfaces/spice.h"
#include "lupnt/interfaces/tai_utc.h"
#include "lupnt/numerics/cheby_fit.h"
#include "lupnt/numerics/math_utils.h"

#define TIME_CONVERSION(from, to, func) \
  {{Time::from, Time::to}, [](Real t) -> Real { return func(t); }}

/// @note
/// D. Folta, N. Bosanac, I. Elliott, L. Mann, R. Mesarch, and J. Rosales,
/// ‘Astrodynamics Convention and Modeling Reference for Lunar, Cislunar, and
/// Libration Point Orbits’, Jan. 2022.
///
/// O. Montenbruck and G. Eberhard, “Satellite Orbits: Models, Methods, and
/// Applications,” Berlin : New York: Springer, 2000.
/// doi: 10.1007/978-3-642-58351-3.

namespace lupnt {

  // Julian Epoch          -4712-01-01T12:00:00 TT
  // Modified Julian Epoch  1858-11-17T00:00:00 TT
  // Fifties Epoch          1950-01-01T00:00:00 TT
  // CCSDS Epoch            1958-01-01T00:00:00 TAI
  // Coordinate Time Epoch  1977-01-01T00:00:00 TAI
  // Galileo Epoch          1999-08-22T00:00:00 UTC
  // GPS Epoch              1980-01-06T00:00:00 UTC
  // J2000 Epoch            2000-01-01T12:00:00 TT

  const std::map<Time, std::string> time_to_string
      = {{Time::UT1, "UT1"}, {Time::UTC, "UTC"}, {Time::TAI, "TAI"}, {Time::TDB, "TDB"},
         {Time::TT, "TT"},   {Time::TCG, "TCG"}, {Time::TCB, "TCB"}, {Time::GPS, "GPS"},
         {Time::TCL, "TCL"}, {Time::LT, "LT"}};

  const std::map<std::string, Time> string_to_time
      = {{"UT1", Time::UT1}, {"UTC", Time::UTC}, {"TAI", Time::TAI}, {"TDB", Time::TDB},
         {"TT", Time::TT},   {"TCG", Time::TCG}, {"TCB", Time::TCB}, {"GPS", Time::GPS},
         {"TCL", Time::TCL}, {"LT", Time::LT}};

  /// @brief Convert time from one time system to another
  /// @param t Time in the original time system
  /// @param from Original time system
  /// @param to Converted time system
  /// @return Real Time in the converted time system
  /// @note
  ///                     TCG
  ///                      |
  /// UT1 -- UTC -- TAI -- TT -> TCB -- TCL -- LT
  ///                |     |      |
  ///               GPS   TDB <---+
  ///
  /// @brief Convert an absolute epoch between time scales.
  ///
  /// Delegates to Epoch, which holds the time-scale relationships.
  ///
  /// PRECISION: the return value is an ABSOLUTE epoch, so it is quantised to the
  /// epoch ULP -- ~245 ns at present-day dates, measured ~38 ns rms / 115 ns peak
  /// for TDB<->TT over 2020-2035. That floor is a property of the return type; no
  /// model improvement can lift it. When the *offset* between two scales is the
  /// quantity of interest, use Epoch (0.0005 ns) or the offset accessors
  /// (TtMinusTdb / TdbMinusTcl / TdbMinusLt).
  Real ConvertTime(Real t, Time from, Time to) {
    if (from == to) return t;
    return Epoch::FromSeconds(t, from).To(to).ToSeconds();
  }

  bool IsCoordinateTimeScale(Time time) {
    switch (time) {
      case Time::TT:
      case Time::TCG:
      case Time::TDB:
      case Time::TCB:
      case Time::TCL:
      case Time::LT: return true;
      default: return false;
    }
  }

  Real ConvertCoordinateTime(Real t, Time from, Time to) {
    LUPNT_CHECK(IsCoordinateTimeScale(from), "Input time scale is not a coordinate time scale",
                "ConvertCoordinateTime");
    LUPNT_CHECK(IsCoordinateTimeScale(to), "Output time scale is not a coordinate time scale",
                "ConvertCoordinateTime");
    return ConvertTime(t, from, to);
  }

  Real ConvertCoordinateTime(Real t, Time from, Time to, const Vec3& x_bcrs) {
    LUPNT_CHECK(IsCoordinateTimeScale(from), "Input time scale is not a coordinate time scale",
                "ConvertCoordinateTime");
    LUPNT_CHECK(IsCoordinateTimeScale(to), "Output time scale is not a coordinate time scale",
                "ConvertCoordinateTime");
    if (from == to) return t;
    if (from == Time::TDB && to == Time::TT) return TDBToTt(t, x_bcrs);
    if (from == Time::TT && to == Time::TDB) return TtToTdb(t, x_bcrs);
    if (from == Time::TT && to == Time::TCB) return TtToTcb(t, x_bcrs);
    if (from == Time::TCB && to == Time::TT) return TcbToTt(t, x_bcrs);
    if (from == Time::TDB && to == Time::TCB) return TtToTcb(TDBToTt(t, x_bcrs), x_bcrs);
    if (from == Time::TCB && to == Time::TDB) return TtToTdb(TcbToTt(t, x_bcrs), x_bcrs);
    if (from == Time::TCG && to == Time::TDB) return TtToTdb(TcgToTt(t), x_bcrs);
    if (from == Time::TCG && to == Time::TCB) return TtToTcb(TcgToTt(t), x_bcrs);
    if (from == Time::TDB && to == Time::TCG) return TtToTcg(TDBToTt(t, x_bcrs));
    if (from == Time::TCB && to == Time::TCG) return TtToTcg(TcbToTt(t, x_bcrs));
    if (from == Time::TCB && to == Time::TCL) return TcbToTcl(t, x_bcrs);
    if (from == Time::TCL && to == Time::TCB) return TclToTcb(t, x_bcrs);
    if (to == Time::TCL) return TcbToTcl(ConvertCoordinateTime(t, from, Time::TCB, x_bcrs), x_bcrs);
    if (from == Time::TCL) return ConvertCoordinateTime(TclToTcb(t, x_bcrs), Time::TCB, to, x_bcrs);
    return ConvertCoordinateTime(t, from, to);
  }

  VecX ConvertCoordinateTime(const VecX& t, Time from, Time to) {
    LUPNT_CHECK(IsCoordinateTimeScale(from), "Input time scale is not a coordinate time scale",
                "ConvertCoordinateTime");
    LUPNT_CHECK(IsCoordinateTimeScale(to), "Output time scale is not a coordinate time scale",
                "ConvertCoordinateTime");
    VecX t_out(t.size());
    for (int i = 0; i < t.size(); ++i) {
      t_out(i) = ConvertCoordinateTime(t(i), from, to);
    }
    return t_out;
  }

  VecX ConvertCoordinateTime(const VecX& t, Time from, Time to, const MatX3& x_bcrs) {
    LUPNT_CHECK(IsCoordinateTimeScale(from), "Input time scale is not a coordinate time scale",
                "ConvertCoordinateTime");
    LUPNT_CHECK(IsCoordinateTimeScale(to), "Output time scale is not a coordinate time scale",
                "ConvertCoordinateTime");
    LUPNT_CHECK(t.size() == x_bcrs.rows(), "Position row count must match time vector size",
                "ConvertCoordinateTime");
    VecX t_out(t.size());
    for (int i = 0; i < t.size(); ++i) {
      t_out(i) = ConvertCoordinateTime(t(i), from, to, x_bcrs.row(i).transpose());
    }
    return t_out;
  }

  VecX ConvertTime(const VecX& t, Time from, Time to) {
    VecX t_out(t.size());
    if (to == Time::TCL) {
      if (from == Time::LT) {
        return LtToTcl(t);
      } else if (from == Time::TCB) {
        return TcbToTcl(t);
      }
      return TcbToTcl(ConvertTime(t, from, Time::TCB));
    }
    if (from == Time::TCL) {
      if (to == Time::LT) {
        return TclToLt(t);
      } else if (to == Time::TCB) {
        return TclToTcb(t);
      }
      return ConvertTime(TclToTcb(t), Time::TCB, to);
    }
    for (int i = 0; i < t.size(); i++) {
      t_out(i) = ConvertTime(t(i), from, to);
    }
    return t_out;
  }

  Real UtcToUt1(Real t_utc) {
    Real mjd_utc = t_utc / SECS_DAY + MJD_J2000_TT;
    Real ut1_utc = GetUt1UtcDifference(mjd_utc);
    Real t_ut1 = t_utc + ut1_utc;
    return t_ut1;
  }

  Real Ut1ToUtc(Real t_ut1) {
    Real mjd_ut1 = t_ut1 / SECS_DAY + MJD_J2000_TT;
    Real ut1_utc = GetUt1UtcDifference(mjd_ut1);
    Real t_utc = t_ut1 - ut1_utc;
    return t_utc;
  }

  Real UtcToTai(Real t_utc) {
    Real mjd_utc = t_utc / SECS_DAY + MJD_J2000_TT;
    Real tai_utc = GetTaiUtcDifference(mjd_utc.val());
    Real t_tai = t_utc + tai_utc;
    return t_tai;
  }

  Real TaiToUtc(Real t_tai) {
    Real mjd_tai = t_tai / SECS_DAY + MJD_J2000_TT;
    Real tai_utc = GetTaiUtcDifference(mjd_tai.val());
    Real t_utc = t_tai - tai_utc;
    return t_utc;
  }

  Real TaiToTt(Real t_tai) { return t_tai + TT_TAI_OFFSET; }

  Real TtToTai(Real t_tt) { return t_tt - TT_TAI_OFFSET; }

  Real TtToTcg(Real t_tt) {
    Real mjd_tt = TimeToMjd(t_tt);
    Real tt_tcg = -L_G / (1.0 - L_G) * (mjd_tt - MJD_COORDINATE_TT_TCG_TCB) * SECS_DAY;
    return t_tt - tt_tcg;
  }

  Real TcgToTt(Real t_tcg) {
    Real mjd_tcg = TimeToMjd(t_tcg);
    Real tt_tcg = -L_G * (mjd_tcg - MJD_COORDINATE_TT_TCG_TCB) * SECS_DAY;
    return t_tcg + tt_tcg;
  }

  namespace {

    constexpr std::array<BodyId, 9> kTtExternalBodies = {
        BodyId::SUN,     BodyId::MERCURY, BodyId::VENUS,  BodyId::MOON,    BodyId::MARS,
        BodyId::JUPITER, BodyId::SATURN,  BodyId::URANUS, BodyId::NEPTUNE,
    };

    Real GetEarthExternalPotentialForTdb(Real t_tdb) {
      Vec3 x_earth = GetBodyPos(t_tdb, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
      Real w_ext = 0.0;
      for (BodyId body : kTtExternalBodies) {
        Vec3 x_body = GetBodyPos(t_tdb, BodyId::SSB, body, Frame::GCRF);
        w_ext += GetBodyGM(body) / (x_earth - x_body).norm();
      }
      return w_ext;
    }

    Real TdbToTtIntegrandC2(Real t_tdb) {
      Vec6 rv_earth = GetBodyPosVel(t_tdb, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
      Real v_earth_2 = rv_earth.tail<3>().squaredNorm();
      Real w_ext = GetEarthExternalPotentialForTdb(t_tdb);
      return 0.5 * v_earth_2 + w_ext;
    }

    Real TdbToTtIntegrandC4(Real t_tdb) {
      Vec6 rv_earth = GetBodyPosVel(t_tdb, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
      Real v_earth_2 = rv_earth.tail<3>().squaredNorm();
      Real w_ext = GetEarthExternalPotentialForTdb(t_tdb);
      return 0.125 * v_earth_2 * v_earth_2 + 1.5 * v_earth_2 * w_ext - 0.5 * w_ext * w_ext;
    }

    std::pair<Real, Real> IntegrateTdbToTtTerms(Real t_tdb) {
      const Real t0_tdb = MjdToTime(MJD_COORDINATE_TT_TCG_TCB) + TDB_0;
      const double t_span = (t_tdb - t0_tdb).val();
      if (std::abs(t_span) < 1.0e-15) return {0.0, 0.0};

      const double sign = t_span >= 0.0 ? 1.0 : -1.0;
      const Real dt_nominal = sign * 0.01 * SECS_DAY;

      Real t_current = t0_tdb;
      Real f2_prev = TdbToTtIntegrandC2(t_current);
      Real f4_prev = TdbToTtIntegrandC4(t_current);
      Real i2 = 0.0;
      Real i4 = 0.0;

      while ((sign > 0.0 && t_current < t_tdb) || (sign < 0.0 && t_current > t_tdb)) {
        Real t_next = t_current + dt_nominal;
        if ((sign > 0.0 && t_next > t_tdb) || (sign < 0.0 && t_next < t_tdb)) {
          t_next = t_tdb;
        }

        Real dt = t_next - t_current;
        Real f2 = TdbToTtIntegrandC2(t_next);
        Real f4 = TdbToTtIntegrandC4(t_next);
        i2 += 0.5 * (f2_prev + f2) * dt;
        i4 += 0.5 * (f4_prev + f4) * dt;

        t_current = t_next;
        f2_prev = f2;
        f4_prev = f4;
      }

      return {i2, i4};
    }

    Real TdbToTtPositionTerm(Real t_tdb, const Vec3& x_bcrs) {
      Vec6 rv_earth = GetBodyPosVel(t_tdb, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
      Vec3 x_earth = rv_earth.head<3>();
      Vec3 v_earth = rv_earth.tail<3>();
      Vec3 r_earth = x_bcrs - x_earth;

      Real v_dot_r = v_earth.dot(r_earth);
      Real v_earth_2 = v_earth.squaredNorm();
      Real w_ext = GetEarthExternalPotentialForTdb(t_tdb);
      const double c2 = C * C;
      const double c4 = c2 * c2;

      return -v_dot_r / c2 - (0.5 * v_earth_2 + 3.0 * w_ext) * v_dot_r / c4;
    }

    // =====================================================================
    // DE440 Eq. (3): TDB - TT relativistic integral
    //
    // R. S. Park, W. M. Folkner, J. G. Williams, D. H. Boggs, "The JPL
    // Planetary and Lunar Ephemerides DE440 and DE441", AJ 161:105 (2021),
    // Eqs. (3)-(7).
    //
    //   TDB - TT = (L_G - L_B)/(1 - L_B) * (TDB - T_0)
    //            + (1 - L_G)/(1 - L_B) * TDB_0
    //            + (1 - L_G)/(1 - L_B) * INT_{T_0+TDB_0}^{TDB} (1/c^2) f2 dt
    //            + (1/c^2) v_E . (r_S - r_E)
    //            - (1 - L_G)/(1 - L_B) * INT_{T_0+TDB_0}^{TDB} (1/c^4) f4 dt
    //            + (1/c^4) (3 w_0E + v_E^2/2) v_E . (r_S - r_E)
    //
    // with   f2 = v_E^2/2 + w_0E + w_LE
    //        f4 = -v_E^4/8 - (3/2) v_E^2 w_0E + 4 v_E . w_iE
    //             + (1/2) w_0E^2 + Delta_E
    //
    // This is the same structure as TdbToTtMinusTdbEq21 above (Turyshev 2026
    // Eq. 21), extended with the three terms that form omits: the solar
    // oblateness potential w_LE (Eq. 5), the mass-weighted velocity potential
    // w_iE (Eq. 6), and Delta_E (Eq. 7). All three are sub-nanosecond over a
    // multi-decade span, so the two forms agree closely; Eq. (3) is provided
    // as the faithful DE440 reference implementation.
    // =====================================================================

    /// Mass parameters for the DE440 Eq. (3) potentials, from Park et al. 2021
    /// Table 2 [m^3/s^2].
    ///
    /// These deliberately use the planetary *system* GM paired with the
    /// system-barycenter position that the DE ephemeris provides for the outer
    /// planets -- that is the potential Earth actually feels from the whole
    /// Jovian/Saturnian/... system.
    ///
    /// They are listed explicitly rather than via GetBodyGM() because LuPNT's
    /// global GM_JUPITER / GM_SATURN / GM_URANUS / GM_NEPTUNE constants are
    /// each exactly 10x smaller than the corresponding *_SYSTEM values (see
    /// constants.h). A planet-only GM is ~0.9996x its system GM -- the moons
    /// are ~4e-4 of the mass -- so those four globals are off by a decimal
    /// place, which would inject a spurious ~ms-level secular drift here.
    struct De440Body {
      BodyId id;
      double gm;
    };
    constexpr std::array<De440Body, 10> kDe440Bodies = {{
        {BodyId::SUN, GM_SUN},
        {BodyId::MERCURY, GM_MERCURY},
        {BodyId::VENUS, GM_VENUS},
        {BodyId::MOON, GM_MOON},
        {BodyId::MARS, GM_MARS_SYSTEM},
        {BodyId::JUPITER, GM_JUPITER_SYSTEM},
        {BodyId::SATURN, GM_SATURN_SYSTEM},
        {BodyId::URANUS, GM_URANUS_SYSTEM},
        {BodyId::NEPTUNE, GM_NEPTUNE_SYSTEM},
        // w_0E sums over ALL bodies other than Earth, so the Pluto system belongs
        // here. Turyshev et al. 2025 Table 4 puts its mean contribution at the EMB
        // at 0.002e-15 = 2e-18 s/s -- ~8% of the residual secular difference
        // between this Eq.(3) implementation and the DE440t TT-TDB ephemeris.
        {BodyId::PLUTO_BARYCENTER, GM_PLUTO_SYSTEM},
    }};

    /// Body states + derived potentials at one epoch, shared by the c^-2 and
    /// c^-4 integrands so the ephemeris is sampled once per step.
    struct De440Terms {
      Vec3 r_E, v_E;
      Real v_E2;
      Real w_0E;  // Eq. (4)  sum_{i!=E} GM_i / r_iE
      Real f2;    // c^-2 integrand
      Real f4;    // c^-4 integrand (as written in Eq. 3)
    };

    /// @brief Newtonian potential [m^2/s^2] at a point in the plane of a uniform
    /// circular ring, at heliocentric distance `r` from the ring centre.
    ///
    ///     U = 2 GM / (pi (a + r)) * K(k),   k = 2 sqrt(a r) / (a + r)
    ///
    /// (Turyshev et al. 2025, ApJ 985:140, Eq. A5 -- the same expression they use
    /// to average a planet's potential at the Earth-Moon barycentre.)
    ///
    /// K is evaluated by the arithmetic-geometric mean, K(k) = pi/(2 AGM(1,k')),
    /// k' = sqrt(1-k^2). AGM converges quadratically (~5 iterations here) and is
    /// differentiable, so this composes with autodiff unlike std::comp_ellint_1.
    Real RingPotential(Real r, double gm_ring, double a_ring) {
      if (gm_ring == 0.0) return Real(0.0);
      Real apr = a_ring + r;
      Real k2 = 4.0 * a_ring * r / (apr * apr);
      Real a = 1.0;
      Real b = sqrt(1.0 - k2);  // complementary modulus k'
      for (int i = 0; i < 12; i++) {
        Real a_next = 0.5 * (a + b);
        b = sqrt(a * b);
        a = a_next;
        if (abs((a - b).val()) < 1e-16) break;
      }
      Real K = PI / (2.0 * a);
      return 2.0 * gm_ring * K / (PI * apr);
    }

    De440Terms ComputeDe440Terms(Real t_tdb) {
      De440Terms out;
      Vec6 rv_E = GetBodyPosVel(t_tdb, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
      out.r_E = rv_E.head<3>();
      out.v_E = rv_E.tail<3>();
      out.v_E2 = out.v_E.squaredNorm();

      // Gather external body states once.
      constexpr int kN = static_cast<int>(kDe440Bodies.size());
      std::array<Vec3, kN> r_i, v_i;
      std::array<Real, kN> gm_i, r_iE;
      for (int i = 0; i < kN; i++) {
        Vec6 rv = GetBodyPosVel(t_tdb, BodyId::SSB, kDe440Bodies[i].id, Frame::GCRF);
        r_i[i] = rv.head<3>();
        v_i[i] = rv.tail<3>();
        gm_i[i] = kDe440Bodies[i].gm;
        r_iE[i] = (out.r_E - r_i[i]).norm();
      }

      // Eq. (4): w_0E = sum_{i!=E} GM_i / r_iE
      out.w_0E = 0.0;
      for (int i = 0; i < kN; i++) out.w_0E += gm_i[i] / r_iE[i];

      // DE440 also integrates 343 asteroids, 30 KBOs and a 36-point Kuiper ring
      // at 44 au. LuPNT has no ephemerides for any of them, so each population is
      // added here as a uniform circular ring, with the masses taken from DE440's
      // OWN integration header (see GM_ASTEROID_BELT / GM_KUIPER_BELT). Those are
      // the values DE440 actually integrated and must NOT be tuned against the
      // DE440t residual.
      {
        const int i_sun_w = 0;  // kDe440Bodies[0].id == BodyId::SUN
        Real r_helio = (out.r_E - r_i[i_sun_w]).norm();
        out.w_0E += RingPotential(r_helio, GM_ASTEROID_BELT, A_ASTEROID_BELT);
        out.w_0E += RingPotential(r_helio, GM_KUIPER_BELT, A_KUIPER_BELT);
      }

      // Eq. (5): solar oblateness potential, with phi_E,sun the heliocentric
      // ecliptic latitude of Earth. Rotate the Sun->Earth vector from the
      // equatorial (ICRF) frame to the ecliptic by the mean J2000 obliquity.
      Real w_LE = 0.0;
      {
        const int i_sun = 0;  // kDe440Bodies[0].id == BodyId::SUN
        Vec3 r_SE = out.r_E - r_i[i_sun];
        Real r_SE_norm = r_SE.norm();
        const double eps = OBLIQUITY_J2000;
        Real z_ecl = -r_SE(1) * std::sin(eps) + r_SE(2) * std::cos(eps);
        Real sin_phi = z_ecl / r_SE_norm;
        w_LE = -GM_SUN * J2_SUN / pow(r_SE_norm, 3) * (R_SUN * R_SUN / 2.0)
               * (3.0 * sin_phi * sin_phi - 1.0);
      }

      // Eq. (6): w_iE = sum_{i!=E} GM_i v_i / r_iE   (vector)
      Vec3 w_iE = Vec3::Zero();
      for (int i = 0; i < kN; i++) w_iE += gm_i[i] / r_iE[i] * v_i[i];

      // Eq. (7): Delta_E. The inner sums run over all bodies (including Earth),
      // and a_i is the Newtonian acceleration of body i (sufficient inside a
      // c^-4 term).
      Real delta_E = 0.0;
      for (int i = 0; i < kN; i++) {
        Real sum_gm_over_r = GM_EARTH / r_iE[i];  // j = Earth
        Vec3 a_i = GM_EARTH * (out.r_E - r_i[i]) / pow(r_iE[i], 3);
        for (int j = 0; j < kN; j++) {
          if (j == i) continue;
          Vec3 d = r_i[j] - r_i[i];
          Real r_ij = d.norm();
          sum_gm_over_r += gm_i[j] / r_ij;
          a_i += gm_i[j] * d / pow(r_ij, 3);
        }
        Vec3 rEi = out.r_E - r_i[i];
        Real proj = v_i[i].dot(rEi) / r_iE[i];
        delta_E += gm_i[i] / r_iE[i]
                   * (-2.0 * v_i[i].squaredNorm() + sum_gm_over_r + 0.5 * proj * proj
                      + 0.5 * a_i.dot(rEi));
      }

      out.f2 = 0.5 * out.v_E2 + out.w_0E + w_LE;
      out.f4 = -0.125 * out.v_E2 * out.v_E2 - 1.5 * out.v_E2 * out.w_0E + 4.0 * out.v_E.dot(w_iE)
               + 0.5 * out.w_0E * out.w_0E + delta_E;
      return out;
    }

    /// Trapezoidal integration of the DE440 Eq. (3) c^-2 and c^-4 integrands
    /// from T_0 + TDB_0 to t_tdb. Returns {INT f2 dt, INT f4 dt}.
    std::pair<Real, Real> IntegrateDe440Terms(Real t_tdb, double step) {
      const Real t0 = MjdToTime(MJD_COORDINATE_TT_TCG_TCB) + TDB_0;
      const double span = (t_tdb - t0).val();
      if (std::abs(span) < 1.0e-15) return {0.0, 0.0};

      const double sign = span >= 0.0 ? 1.0 : -1.0;
      const Real dt_nominal = sign * step;

      Real t_current = t0;
      De440Terms prev = ComputeDe440Terms(t_current);
      Real i2 = 0.0, i4 = 0.0;

      while ((sign > 0.0 && t_current < t_tdb) || (sign < 0.0 && t_current > t_tdb)) {
        Real t_next = t_current + dt_nominal;
        if ((sign > 0.0 && t_next > t_tdb) || (sign < 0.0 && t_next < t_tdb)) t_next = t_tdb;
        Real dt = t_next - t_current;
        De440Terms cur = ComputeDe440Terms(t_next);
        i2 += 0.5 * (prev.f2 + cur.f2) * dt;
        i4 += 0.5 * (prev.f4 + cur.f4) * dt;
        t_current = t_next;
        prev = cur;
      }
      return {i2, i4};
    }

    Real TdbToTtMinusTdbEq21(Real t_tdb, const Vec3& x_bcrs) {
      auto [i2, i4] = IntegrateTdbToTtTerms(t_tdb);
      const Real t0 = MjdToTime(MJD_COORDINATE_TT_TCG_TCB);
      const double c2 = C * C;
      const double c4 = c2 * c2;
      const Real position = -TdbToTtPositionTerm(t_tdb, x_bcrs);
      const Real bracket = TDB_0 + i2 / c2 + i4 / c4 + position;
      return (L_B - L_G) / (1.0 - L_B) * (t_tdb - t0) - (1.0 - L_G) / (1.0 - L_B) * bracket;
    }

  }  // namespace

  // Optional piecewise-Chebyshev fit of TT-TDB [s] as a function of TDB,
  // sampled from the high-fidelity DE440t TT-TDB ephemeris (see
  // InitTtMinusTdbFit). When populated, TDBToTt/TtToTdb use it for fast,
  // SPICE-free evaluation (important for OpenMP-parallel code, since the
  // SPICE interface serializes on an omp-critical section); otherwise they
  // fall back to the analytic series below.
  namespace {
    std::optional<ChebyshevFitModel> s_tt_minus_tdb_fit;

    // Evaluate the fitted TT-TDB [s] at epoch `t` (TDB, or TT to within the
    // ~2 ms |TT-TDB| offset). Returns false if no fit covers `t`.
    bool EvalTtMinusTdbFit(Real t, Real* tt_minus_tdb) {
      VecX val(1);
      bool ok = false;
#pragma omp critical
      {
        ok = s_tt_minus_tdb_fit.has_value() && s_tt_minus_tdb_fit->Eval(t, &val, nullptr);
      }
      if (ok) *tt_minus_tdb = val(0);
      return ok;
    }

    // Selects which model TDBToTt/TtToTdb use when no Chebyshev fit covers the
    // epoch and auto-fitting is unavailable.
    TtTdbModel s_tt_tdb_model = TtTdbModel::FITTED;
    double s_de440_step = 0.01 * SECS_DAY;  // integration step for Eq. (3) [s]

    // Auto-fit: build a DE440t Chebyshev fit on demand so the DEFAULT accuracy
    // is ~0.4 ps rather than the ~17 us of the analytic series. Disable with
    // SetTtTdbAutoFit(false) to get the historical behaviour.
    bool s_tt_tdb_autofit = true;
    bool s_tt_tdb_autofit_unavailable = false;  // sticky: kernel missing, stop retrying

    /// Fit a decade-wide window covering `t`, if one is not already loaded.
    ///
    /// The window is snapped to a fixed decade grid so that sequential epochs
    /// reuse a single fit instead of refitting at every step. Returns true if a
    /// fit covering `t` should now exist. Never throws: if the DE440t TT-TDB
    /// segment is unavailable the failure is recorded and the caller falls back
    /// to the analytic series.
    ///
    /// MUST NOT be called from inside an omp critical region -- it samples SPICE
    /// (which takes the unnamed critical lock) and then stores the fit.
    bool TryAutoFitTtMinusTdb(Real t) {
      bool enabled = false;
#pragma omp critical
      {
        enabled = s_tt_tdb_autofit && !s_tt_tdb_autofit_unavailable;
      }
      if (!enabled) return false;

      constexpr double kYear = 365.25 * SECS_DAY;
      constexpr double kWindow = 10.0 * kYear;
      const double centre = std::floor(t.val() / kWindow) * kWindow;
      const double t0 = centre - kYear;  // 1 yr of margin either side
      const double t1 = centre + kWindow + kYear;

      try {
        InitTtMinusTdbFit(Real(t0), Real(t1));
      } catch (...) {
#pragma omp critical
        {
          s_tt_tdb_autofit_unavailable = true;
        }
        return false;
      }
      return true;
    }
  }  // namespace

  void SetTtTdbModel(TtTdbModel model) {
#pragma omp critical
    {
      s_tt_tdb_model = model;
    }
  }

  TtTdbModel GetTtTdbModel() {
    TtTdbModel m;
#pragma omp critical
    {
      m = s_tt_tdb_model;
    }
    return m;
  }

  void SetTtTdbAutoFit(bool enable) {
#pragma omp critical
    {
      s_tt_tdb_autofit = enable;
      if (enable) s_tt_tdb_autofit_unavailable = false;  // allow a retry
    }
  }

  bool GetTtTdbAutoFit() {
    bool v;
#pragma omp critical
    {
      v = s_tt_tdb_autofit;
    }
    return v;
  }

  void SetDe440TtTdbStep(double step_s) {
    LUPNT_CHECK(step_s > 0.0, "DE440 TDB-TT integration step must be positive", "TimeConversions");
#pragma omp critical
    {
      s_de440_step = step_s;
    }
  }

  double GetDe440TtTdbStep() {
    double s;
#pragma omp critical
    {
      s = s_de440_step;
    }
    return s;
  }

  /// @brief TDB - TT [s] at BCRS station position `x_bcrs`, from DE440 Eq. (3).
  ///
  /// Park et al. 2021 (AJ 161:105), Eqs. (3)-(7). Integrates the c^-2 and c^-4
  /// integrands from T_0 + TDB_0 to `t_tdb` with a trapezoidal rule.
  Real TdbMinusTtDe440(Real t_tdb, const Vec3& x_bcrs) {
    const double step = GetDe440TtTdbStep();
    auto [i2, i4] = IntegrateDe440Terms(t_tdb, step);
    const Real t0 = MjdToTime(MJD_COORDINATE_TT_TCG_TCB);
    const double c2 = C * C;
    const double c4 = c2 * c2;
    const double sB = (1.0 - L_G) / (1.0 - L_B);

    // Endpoint (station) terms: v_E . (r_S - r_E) evaluated at t_tdb.
    De440Terms end = ComputeDe440Terms(t_tdb);
    Vec3 r_rel = x_bcrs - end.r_E;
    Real v_dot_r = end.v_E.dot(r_rel);

    return (L_G - L_B) / (1.0 - L_B) * (t_tdb - t0) + sB * TDB_0 + sB * i2 / c2 + v_dot_r / c2
           - sB * i4 / c4 + (3.0 * end.w_0E + 0.5 * end.v_E2) * v_dot_r / c4;
  }

  /// @brief TDB - TT [s] at the geocenter (r_S = r_E), from DE440 Eq. (3).
  ///
  /// At the geocenter r_S = r_E, so both endpoint (station) terms vanish and
  /// only the two integrals remain.
  Real TdbMinusTtDe440(Real t_tdb) {
    const double step = GetDe440TtTdbStep();
    auto [i2, i4] = IntegrateDe440Terms(t_tdb, step);
    const Real t0 = MjdToTime(MJD_COORDINATE_TT_TCG_TCB);
    const double c2 = C * C;
    const double c4 = c2 * c2;
    const double sB = (1.0 - L_G) / (1.0 - L_B);
    return (L_G - L_B) / (1.0 - L_B) * (t_tdb - t0) + sB * TDB_0 + sB * i2 / c2 - sB * i4 / c4;
  }

  /// @brief Vectorized geocentric DE440 Eq. (3) TDB - TT [s].
  ///
  /// Evaluating the scalar overload at N epochs would re-integrate from
  /// T_0 + TDB_0 N times. This instead sorts the epochs and performs a single
  /// forward sweep (plus one backward sweep for epochs before T_0), recording
  /// the accumulated integrals as each epoch is passed -- the same strategy
  /// used by InitLtMinusTtFit.
  VecX TdbMinusTtDe440(const VecX& t_tdb) {
    const int n = static_cast<int>(t_tdb.size());
    VecX out(n);
    if (n == 0) return out;

    const double step = GetDe440TtTdbStep();
    const double t0 = (MjdToTime(MJD_COORDINATE_TT_TCG_TCB) + TDB_0).val();
    const Real t_ref = MjdToTime(MJD_COORDINATE_TT_TCG_TCB);
    const double c2 = C * C;
    const double c4 = c2 * c2;
    const double sB = (1.0 - L_G) / (1.0 - L_B);

    auto finish = [&](double t, double i2, double i4) {
      return (L_G - L_B) / (1.0 - L_B) * (t - t_ref.val()) + sB * TDB_0 + sB * i2 / c2
             - sB * i4 / c4;
    };

    std::vector<int> order(n);
    for (int i = 0; i < n; i++) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return t_tdb(a).val() < t_tdb(b).val(); });

    // Forward sweep over epochs at/after T_0 + TDB_0.
    {
      double i2 = 0.0, i4 = 0.0, t_cur = t0;
      De440Terms prev = ComputeDe440Terms(t_cur);
      for (int k : order) {
        double tk = t_tdb(k).val();
        if (tk < t0) continue;
        while (t_cur < tk) {
          double t_next = std::min(t_cur + step, tk);
          double dt = t_next - t_cur;
          De440Terms cur = ComputeDe440Terms(t_next);
          i2 += (0.5 * (prev.f2 + cur.f2) * dt).val();
          i4 += (0.5 * (prev.f4 + cur.f4) * dt).val();
          t_cur = t_next;
          prev = cur;
        }
        out(k) = finish(tk, i2, i4);
      }
    }

    // Backward sweep over epochs before T_0 + TDB_0 (descending).
    {
      double i2 = 0.0, i4 = 0.0, t_cur = t0;
      De440Terms prev = ComputeDe440Terms(t_cur);
      for (auto it = order.rbegin(); it != order.rend(); ++it) {
        double tk = t_tdb(*it).val();
        if (tk >= t0) continue;
        while (t_cur > tk) {
          double t_next = std::max(t_cur - step, tk);
          double dt = t_next - t_cur;  // negative
          De440Terms cur = ComputeDe440Terms(t_next);
          i2 += (0.5 * (prev.f2 + cur.f2) * dt).val();
          i4 += (0.5 * (prev.f4 + cur.f4) * dt).val();
          t_cur = t_next;
          prev = cur;
        }
        out(*it) = finish(tk, i2, i4);
      }
    }
    return out;
  }

  /// @brief
  /// @param t_tdb
  /// @return
  /// @ref
  /// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/time.html#The%20Relationship%20between%20TT%20and%20TDB
  namespace {
    // The two-term IAU / IERS TN36 series (as used by Orekit and GMAT):
    //   TDB - TT = 0.001658 sin(g) + 0.000014 sin(2g),  g = Earth mean anomaly.
    Real TtMinusTdbSeries(Real t_tdb) {
      Real days_j2000 = t_tdb / SECS_DAY;
      Real g = (357.53 + 0.9856003 * days_j2000) * RAD;
      return -(0.001658 * sin(g) + 0.000014 * sin(2.0 * g));
    }
  }  // namespace

  /// @brief TT - TDB [s] under the active model, returned as a *small offset*.
  ///
  /// Prefer this over `TDBToTt(t) - t` whenever the offset itself is the
  /// quantity of interest. TT-TDB is ~1e-4 s while an absolute epoch is ~1e9 s
  /// at present-day dates, so forming `t + offset` and subtracting `t` back off
  /// snaps the result to the epoch's ULP grid (~0.25 us at 2030). Returning the
  /// offset directly keeps full double precision (~1e-20 s).
  Real TtMinusTdb(Real t_tdb) {
    // The selected model decides which computation runs, so an explicit choice
    // is never silently overridden by a cached fit.
    switch (GetTtTdbModel()) {
      case TtTdbModel::ANALYTIC: return TtMinusTdbSeries(t_tdb);
      case TtTdbModel::DE440_INTEGRAL:
        // DE440 Eq. (3) integral -- integrates from T_0, so much slower.
        return -TdbMinusTtDe440(t_tdb);
      case TtTdbModel::FITTED: break;  // handled below
    }

    // FITTED (default): a fit covering this epoch, an on-demand fit if auto-fit
    // is enabled, or the series as a last resort when no fit can be produced.
    Real tt_minus_tdb;
    if (EvalTtMinusTdbFit(t_tdb, &tt_minus_tdb)) return tt_minus_tdb;
    if (TryAutoFitTtMinusTdb(t_tdb) && EvalTtMinusTdbFit(t_tdb, &tt_minus_tdb)) return tt_minus_tdb;
    return TtMinusTdbSeries(t_tdb);
  }

  Real TDBToTt(Real t_tdb) { return t_tdb + TtMinusTdb(t_tdb); }

  /// @brief Convert TDB to TT for an event at BCRS position x_bcrs.
  ///
  /// Uses the Turyshev 2026 Eq. (21)/(22) integral form, including the
  /// position-dependent Earth term where r_E = x - x_E and x is the event
  /// position in BCRS.
  Real TDBToTt(Real t_tdb, const Vec3& x_bcrs) {
    return t_tdb + TdbToTtMinusTdbEq21(t_tdb, x_bcrs);
  }

  /// @brief
  /// @param t_tt
  /// @return
  /// @ref
  /// Astrodynamics Convention and Modeling Reference for Lunar, Cislunar, and Libration Point
  /// Orbits
  /// @brief TT -> TDB.
  ///
  /// Delegates to TtMinusTdb() so both directions share one model chain (fit ->
  /// DE440 integral -> auto-fit -> analytic series); this is what makes
  /// TDBToTt(TtToTdb(t)) == t exact regardless of which model is active.
  ///
  /// The fit and the integral are parameterized by TDB, but |TT-TDB| < 2 ms with
  /// rate ~3e-10, so evaluating at t_tt costs < 1e-12 s.
  Real TtToTdb(Real t_tt) { return t_tt - TtMinusTdb(t_tt); }

  Real TtToTdb(Real t_tt, const Vec3& x_bcrs) {
    Real t_tdb = TtToTdb(t_tt);
    for (int iter = 0; iter < 10; ++iter) {
      Real delta = t_tt - TDBToTt(t_tdb, x_bcrs);
      t_tdb += delta;
      if (std::abs(delta.val()) < 1.0e-12) break;
    }
    return t_tdb;
  }

  Real TaiToGps(Real t_tai) {
    Real t_gps = t_tai - 19.0;
    return t_gps;
  }

  Real GpsToTai(Real t_gps) {
    Real t_tai = t_gps + 19.0;
    return t_tai;
  }

  Real TcbToTdb(Real t_tcb) {
    Real mjd_tcb = TimeToMjd(t_tcb);
    Real t_tdb = t_tcb - L_B * (mjd_tcb - MJD_COORDINATE_TT_TCG_TCB) * SECS_DAY + TDB_0;
    return t_tdb;
  }

  Real TtToTcb(Real t_tt) {
    Real t_tdb = TtToTdb(t_tt);
    const Real t0_tcb = MjdToTime(MJD_COORDINATE_TT_TCG_TCB);
    return (t_tdb - L_B * t0_tcb - TDB_0) / (1.0 - L_B);
  }

  Real TtToTcb(Real t_tt, const Vec3& x_bcrs) {
    Real t_tdb = TtToTdb(t_tt, x_bcrs);
    const Real t0_tcb = MjdToTime(MJD_COORDINATE_TT_TCG_TCB);
    return (t_tdb - L_B * t0_tcb - TDB_0) / (1.0 - L_B);
  }

  Real TcbToTt(Real t_tcb, const Vec3& x_bcrs) { return TDBToTt(TcbToTdb(t_tcb), x_bcrs); }

  Real EarthRotationAngle(Real t_ut1) {
    double theta_0 = 0.7790572732640;
    double dtheta_dt = 1.00273781191135448;
    Real theta_era = TWO_PI * (theta_0 + dtheta_dt * t_ut1 / SECS_DAY);
    return WrapToPi(theta_era);
  }

  Real GregorianToMjd(int year, int month, int day, int hour, int min, Real sec) {
    if (month <= 2) {
      month += 12;
      year--;
    }
    int b;
    if ((10000L * year + 100L * month + day) <= 15821004L)
      b = -2 + ((year + 4716) / 4) - 1179;  // Julian calendar
    else
      b = (year / 400) - (year / 100) + (year / 4);  // Gregorian calendar

    double mjd_midnight = 365L * year - 679004L + b + int(30.6001 * (month + 1)) + day;
    Real frac_of_day = (hour + min / 60.0 + sec / 3600.0) / 24.0;
    return mjd_midnight + frac_of_day;
  }

  std::tuple<int, int, int, int, int, Real> MjdToGregorian(Real mjd) {
    long a, b, c, d, e, f;
    a = long(mjd + 2400001.0);  // Convert Julian day number to calendar date
    if (a < 2299161) {          // Julian calendar
      b = 0;
      c = a + 1524;
    } else {  // Gregorian calendar
      b = long((a - 1867216.25) / 36524.25);
      c = a + b - (b / 4) + 1525;
    }
    d = long((c - 122.1) / 365.25);
    e = 365 * d + d / 4;
    f = long((c - e) / 30.6001);
    int day = c - e - int(30.6001 * f);
    int month = f - 1 - 12 * (f / 14);
    int year = d - 4715 - ((7 + month) / 10);
    Real hours = HOURS_DAY * (mjd - floor(mjd));
    int hour = int(hours);
    Real x = (hours - hour) * MINS_HOUR;
    int min = int(x);
    Real sec = (x - min) * SECS_MINUTE;
    return std::make_tuple(year, month, day, hour, min, sec);
  }

  Real GregorianToTime(int year, int month, int day, int hour, int min, Real sec) {
    Real mjd = GregorianToMjd(year, month, day, hour, min, sec);
    return MjdToTime(mjd);
  }

  Real GregorianToTime(const std::string& date) {
    int year, month, day, hour, min;
    double sec;
    sscanf(date.c_str(), "%d-%d-%dT%d:%d:%lf", &year, &month, &day, &hour, &min, &sec);
    return GregorianToTime(year, month, day, hour, min, sec);
  }

  /// @brief Greenwich Mean Sidereal Time
  /// @param mjd_ut1 UT1 (Modified Julian Date)
  /// @return GMST [rad]
  Real GreenwichMeanSiderealTime(Real mjd_ut1) {
    Real mjd0 = floor(mjd_ut1);
    Real ut1 = SECS_DAY * (mjd_ut1 - mjd0);  // [s]
    Real T0 = (mjd0 - MJD_J2000_TT) / DAYS_CENTURY;
    Real T = (mjd_ut1 - MJD_J2000_TT) / DAYS_CENTURY;

    Real gmst = 24110.54841 + 8640184.812866 * T0 + 1.002737909350795 * ut1
                + (0.093104 - 6.2e-6 * T) * T * T;  // [s]

    return TWO_PI * frac(gmst / SECS_DAY);  // [rad]
  }

  Real MjdToTime(Real mjd) { return (mjd - MJD_J2000_TT) * SECS_DAY; }

  Real TimeToMjd(Real t) { return t / SECS_DAY + MJD_J2000_TT; }

  Real JdToTime(Real jd) { return (jd - JD_J2000_TT) * SECS_DAY; }

  Real TimeToJd(Real t) { return t / SECS_DAY + JD_J2000_TT; }

  /// @brief Convert Modified Julian Date to date string
  /// @param mjd Modified Julian Date
  /// @param precision Number of seconds precision
  /// @return Date string
  std::string MjdToGregorianString(Real mjd, int precision) {
    double pow10 = pow(10, precision);
    Real mjd_round = (round(mjd * SECS_DAY * pow10, precision) + 0.1) / (SECS_DAY * pow10);
    auto [year, month, day, hour, min, sec] = MjdToGregorian(mjd_round);
    std::stringstream ss;
    sec = round(sec, precision);
    ss << year << "-";
    ss << std::setw(2) << std::setfill('0') << month << "-";
    ss << std::setw(2) << std::setfill('0') << day << "T";
    ss << std::setw(2) << std::setfill('0') << hour << ":";
    ss << std::setw(2) << std::setfill('0') << min << ":";
    ss << std::setw(2) << std::setfill('0') << floor(sec);
    if (precision > 0) {
      ss << "." << std::fixed << std::setprecision(0) << std::setw(precision) << std::setfill('0')
         << floor((sec - floor(sec)) * pow(10, precision));
    }
    return ss.str();
  }

  std::string TimeToGregorianString(Real t, int precision) {
    Real mjd = TimeToMjd(t);
    return MjdToGregorianString(mjd, precision);
  }

  /// @brief Greenwich Apparent Sidereal Time
  /// @param mjd_ut1 UT1 (Modified Julian Date)
  /// @return GAST [rad]
  Real GreenwichApparentSiderealTime(Real mjd_ut1) {
    return mod(GreenwichMeanSiderealTime(mjd_ut1) + EquinoxEquation(mjd_ut1), TWO_PI);
  }

  namespace {

    // Bodies (other than the Moon) contributing to the selenocentric external
    // potential, with their DE440 mass parameters [m^3/s^2].
    //
    // These use the planetary *system* GM paired with the system-barycenter
    // position the DE ephemeris returns for the outer planets -- that is the
    // potential the Moon actually feels from the whole Jovian/Saturnian/...
    // system, and it matches what the reference lunar time ephemerides assume.
    // Using the planet-only GM here instead leaves the moons' mass out, which
    // accumulates as a ~15 ns/yr drift in TDB-TCL (validated against
    // LTE_DE440_TDBmTCL.bsp).
    constexpr std::array<De440Body, 10> kTclExternalBodies = {{
        {BodyId::SUN, GM_SUN},
        {BodyId::MERCURY, GM_MERCURY},
        {BodyId::VENUS, GM_VENUS},
        {BodyId::EARTH, GM_EARTH},
        {BodyId::MARS, GM_MARS_SYSTEM},
        {BodyId::JUPITER, GM_JUPITER_SYSTEM},
        {BodyId::SATURN, GM_SATURN_SYSTEM},
        {BodyId::URANUS, GM_URANUS_SYSTEM},
        {BodyId::NEPTUNE, GM_NEPTUNE_SYSTEM},
        // Kept in step with kDe440Bodies: the external potential sums over all
        // bodies other than the Moon, so the Pluto system belongs here too.
        {BodyId::PLUTO_BARYCENTER, GM_PLUTO_SYSTEM},
    }};

    Vec6 GetMoonBcrsStateForTcb(Real t_tcb) {
      return GetBodyPosVel(TcbToTdb(t_tcb), BodyId::SSB, BodyId::MOON, Frame::GCRF);
    }

    Real GetLunarExternalPotentialForTcb(Real t_tcb) {
      Real t_tdb = TcbToTdb(t_tcb);
      Vec3 x_moon = GetBodyPos(t_tdb, BodyId::SSB, BodyId::MOON, Frame::GCRF);
      Real w_ext = 0.0;

      Vec3 x_sun_bcrs = GetBodyPos(t_tdb, BodyId::SSB, BodyId::SUN, Frame::GCRF);
      for (const De440Body& b : kTclExternalBodies) {
        Vec3 x_body = GetBodyPos(t_tdb, BodyId::SSB, b.id, Frame::GCRF);
        w_ext += b.gm / (x_moon - x_body).norm();
      }
      // Same small-body ring populations as w_0E in the Eq.(3) path, evaluated at
      // the Moon's heliocentric distance. Masses from DE440's header; see
      // constants.h. Unlike the Earth-Moon difference form used for LT-TT, this
      // is an ABSOLUTE potential, so the distant populations do not cancel and
      // genuinely belong here.
      {
        Real r_helio_moon = (x_moon - x_sun_bcrs).norm();
        w_ext += RingPotential(r_helio_moon, GM_ASTEROID_BELT, A_ASTEROID_BELT);
        w_ext += RingPotential(r_helio_moon, GM_KUIPER_BELT, A_KUIPER_BELT);
      }
      return w_ext;
    }

    Real TcbToTclIntegrandC2(Real t_tcb) {
      Vec6 rv_moon = GetMoonBcrsStateForTcb(t_tcb);
      Real v_moon_2 = rv_moon.tail<3>().squaredNorm();
      Real w_ext = GetLunarExternalPotentialForTcb(t_tcb);
      return 0.5 * v_moon_2 + w_ext;
    }

    Real TcbToTclIntegrandC4(Real t_tcb) {
      Vec6 rv_moon = GetMoonBcrsStateForTcb(t_tcb);
      Real v_moon_2 = rv_moon.tail<3>().squaredNorm();
      Real w_ext = GetLunarExternalPotentialForTcb(t_tcb);
      return 0.125 * v_moon_2 * v_moon_2 + 1.5 * v_moon_2 * w_ext - 0.5 * w_ext * w_ext;
    }

    std::pair<Real, Real> IntegrateTcbToTclTerms(Real t_tcb) {
      const Real t0_tcb = MjdToTime(MJD_COORDINATE_TT_TCG_TCB);
      const double t_span = (t_tcb - t0_tcb).val();
      if (std::abs(t_span) < 1.0e-15) return {0.0, 0.0};

      const double sign = t_span >= 0.0 ? 1.0 : -1.0;
      const Real dt_nominal = sign * 0.01 * SECS_DAY;

      Real t_current = t0_tcb;
      Real f2_prev = TcbToTclIntegrandC2(t_current);
      Real f4_prev = TcbToTclIntegrandC4(t_current);
      Real i2 = 0.0;
      Real i4 = 0.0;

      while ((sign > 0.0 && t_current < t_tcb) || (sign < 0.0 && t_current > t_tcb)) {
        Real t_next = t_current + dt_nominal;
        if ((sign > 0.0 && t_next > t_tcb) || (sign < 0.0 && t_next < t_tcb)) {
          t_next = t_tcb;
        }

        Real dt = t_next - t_current;
        Real f2 = TcbToTclIntegrandC2(t_next);
        Real f4 = TcbToTclIntegrandC4(t_next);
        i2 += 0.5 * (f2_prev + f2) * dt;
        i4 += 0.5 * (f4_prev + f4) * dt;

        t_current = t_next;
        f2_prev = f2;
        f4_prev = f4;
      }

      return {i2, i4};
    }

    Real TcbToTclPositionTerm(Real t_tcb, const Vec3& x_bcrs) {
      Vec6 rv_moon = GetMoonBcrsStateForTcb(t_tcb);
      Vec3 x_moon = rv_moon.head<3>();
      Vec3 v_moon = rv_moon.tail<3>();
      Vec3 r_moon = x_bcrs - x_moon;

      Real v_dot_r = v_moon.dot(r_moon);
      Real v_moon_2 = v_moon.squaredNorm();
      Real w_ext = GetLunarExternalPotentialForTcb(t_tcb);
      const double c2 = C * C;
      const double c4 = c2 * c2;

      return -v_dot_r / c2 - (0.5 * v_moon_2 + 3.0 * w_ext) * v_dot_r / c4;
    }

  }  // namespace

  /// @brief Convert Barycentric Coordinate Time (TCB) to Lunar Coordinate Time (TCL).
  ///
  /// Implements Turyshev 2026, Eq. (23)/(25), using the lunar center as the BCRS
  /// clock position. With x = x_M, the position term v_M · r_M is zero.
  Real TcbToTcl(Real t_tcb) {
    Vec3 x_moon = GetMoonBcrsStateForTcb(t_tcb).head<3>();
    return TcbToTcl(t_tcb, x_moon);
  }

  /// @brief Convert TCB to TCL for a clock at BCRS position x_bcrs.
  ///
  /// Turyshev Eq. (23)/(25) uses r_M = x - x_M, where x is the clock position in
  /// BCRS, x_M is the Moon's BCRS position, v_M = dx_M/dTCB, and
  /// w_ext(x_M) = sum_{B != M} GM_B / |x_M - x_B|.
  Real TcbToTcl(Real t_tcb, const Vec3& x_bcrs) {
    auto [i2, i4] = IntegrateTcbToTclTerms(t_tcb);
    const double c2 = C * C;
    const double c4 = c2 * c2;
    return t_tcb - i2 / c2 - i4 / c4 + TcbToTclPositionTerm(t_tcb, x_bcrs);
  }

  /// @brief TDB - TCL [s] at the lunar centre, returned as a *small offset*.
  ///
  /// Prefer this over `t - ConvertTime(t, TDB, TCL)`. Although TDB-TCL is ~1 s
  /// (so the epoch ULP is only a ~2e-7 relative effect), differencing two ~1e9 s
  /// epochs still snaps the result to a ~0.25 us grid -- far above the
  /// nanosecond level at which this conversion can be validated. Assembling the
  /// offset from its parts keeps full precision:
  ///
  ///   TDB - TCL = (TDB - TCB) + (TCB - TCL)
  ///             = [-L_B (TCB - T_0) + TDB_0] + [i2/c^2 + i4/c^4]
  ///
  /// The first bracket multiplies a large elapsed time by L_B ~ 1.6e-8, so its
  /// own ULP contributes only ~6e-15 s. The position term vanishes at the lunar
  /// centre.
  namespace {
    std::optional<ChebyshevFitModel> s_tdb_minus_tcl_fit;

    bool EvalTdbMinusTclFit(Real t_tdb, Real* out) {
      VecX val(1);
      bool ok = false;
#pragma omp critical
      {
        ok = s_tdb_minus_tcl_fit.has_value() && s_tdb_minus_tcl_fit->Eval(t_tdb, &val, nullptr);
      }
      if (ok) *out = val(0);
      return ok;
    }

    bool s_tdb_tcl_autofit = true;
    bool s_tdb_tcl_autofit_unavailable = false;

    /// Fit a decade-wide window covering `t_tdb`, snapped to a fixed grid so
    /// sequential epochs reuse one fit. Never throws. MUST NOT be called from
    /// inside an omp critical region.
    bool TryAutoFitTdbMinusTcl(Real t_tdb) {
      bool enabled = false;
#pragma omp critical
      {
        enabled = s_tdb_tcl_autofit && !s_tdb_tcl_autofit_unavailable;
      }
      if (!enabled) return false;

      constexpr double kYear = 365.25 * SECS_DAY;
      constexpr double kWindow = 10.0 * kYear;
      const double t0_tdb = MjdToTime(MJD_COORDINATE_TT_TCG_TCB).val();
      // The integral is anchored at T_0; epochs before it must use the slow path.
      if (t_tdb.val() < t0_tdb) return false;
      const double centre = std::floor(t_tdb.val() / kWindow) * kWindow;
      const double win_start = std::max(centre - kYear, t0_tdb);
      try {
        // 4-day segments: TDB-TCL carries the 27.32-day lunar month, which
        // 16-day segments resolve poorly (0.17 ns interpolation error).
        InitTdbMinusTclFit(Real(win_start), Real(centre + kWindow + kYear), 4.0 * SECS_DAY, 13);
      } catch (...) {
#pragma omp critical
        {
          s_tdb_tcl_autofit_unavailable = true;
        }
        return false;
      }
      return true;
    }
  }  // namespace

  Real TdbMinusTcl(Real t_tdb) {
    // Fast path: fitted model. The slow path below integrates from T_0 (1977) on
    // every call (~12 s), so without a fit any per-step lunar conversion is
    // unusable. Auto-fitting builds one decade-wide window on demand.
    Real fitted;
    if (EvalTdbMinusTclFit(t_tdb, &fitted)) return fitted;
    if (TryAutoFitTdbMinusTcl(t_tdb) && EvalTdbMinusTclFit(t_tdb, &fitted)) return fitted;

    const Real t0 = MjdToTime(MJD_COORDINATE_TT_TCG_TCB);
    // TDB -> TCB (inverse of TcbToTdb).
    Real t_tcb = (t_tdb - L_B * t0 - TDB_0) / (1.0 - L_B);
    auto [i2, i4] = IntegrateTcbToTclTerms(t_tcb);
    const double c2 = C * C;
    const double c4 = c2 * c2;
    return -L_B * (t_tcb - t0) + TDB_0 + i2 / c2 + i4 / c4;
  }

  /// @brief TDB - TL [s], returned as a *small offset*.
  ///
  /// Assembled from two small quantities, so no absolute epoch is ever formed:
  ///   TL - TDB = (TL - TT) + (TT - TDB)
  /// Uses whichever TT-TDB model is active (fit / DE440 integral / series).
  ///
  /// This is one of two independent routes to LT. The other, TDB -> TCL -> LT
  /// through Epoch's conversion graph, reaches the same instant via the absolute
  /// lunar potential (TdbMinusTcl, which sums every body plus the small-body
  /// rings) rather than the Earth-Moon difference form of Eq. (73), in which
  /// distant bodies cancel and never appear.
  ///
  /// The two agree at T_0 and diverge by ~0.06 ns/yr thereafter: 0.004 ns half a
  /// year after T_0, 3.5 ns by 2035. Fitting the difference against (t - T_0)
  /// gives an intercept of 0.01 ns and a residual of 0.06 ns, i.e. it is a pure
  /// accumulation of a small constant difference between the two integrands --
  /// not an offset and not an interpolation artefact. The rate is the same order
  /// as the residual TT-TDB drift against DE440t (-0.046 ns/yr).
  Real TdbMinusLt(Real t_tdb) { return -(TdbToLtMinusTt(t_tdb) + TtMinusTdb(t_tdb)); }

  VecX TcbToTcl(const VecX& t_tcb) {
    VecX t_tcl(t_tcb.size());
    for (int i = 0; i < t_tcb.size(); ++i) {
      t_tcl(i) = TcbToTcl(t_tcb(i));
    }
    return t_tcl;
  }

  VecX TcbToTcl(const VecX& t_tcb, const MatX3& x_bcrs) {
    LUPNT_CHECK(t_tcb.size() == x_bcrs.rows(), "Position row count must match time vector size",
                "TcbToTcl");
    VecX t_tcl(t_tcb.size());
    for (int i = 0; i < t_tcb.size(); ++i) {
      t_tcl(i) = TcbToTcl(t_tcb(i), x_bcrs.row(i).transpose());
    }
    return t_tcl;
  }

  Real TclToTcb(Real t_tcl) {
    Real t_tcb = t_tcl;
    for (int iter = 0; iter < 10; ++iter) {
      Real delta = t_tcl - TcbToTcl(t_tcb);
      t_tcb += delta;
      if (std::abs(delta.val()) < 1.0e-12) break;
    }
    return t_tcb;
  }

  Real TclToTcb(Real t_tcl, const Vec3& x_bcrs) {
    Real t_tcb = t_tcl;
    for (int iter = 0; iter < 10; ++iter) {
      Real delta = t_tcl - TcbToTcl(t_tcb, x_bcrs);
      t_tcb += delta;
      if (std::abs(delta.val()) < 1.0e-12) break;
    }
    return t_tcb;
  }

  VecX TclToTcb(const VecX& t_tcl) {
    VecX t_tcb(t_tcl.size());
    for (int i = 0; i < t_tcl.size(); ++i) {
      t_tcb(i) = TclToTcb(t_tcl(i));
    }
    return t_tcb;
  }

  VecX TclToTcb(const VecX& t_tcl, const MatX3& x_bcrs) {
    LUPNT_CHECK(t_tcl.size() == x_bcrs.rows(), "Position row count must match time vector size",
                "TclToTcb");
    VecX t_tcb(t_tcl.size());
    for (int i = 0; i < t_tcl.size(); ++i) {
      t_tcb(i) = TclToTcb(t_tcl(i), x_bcrs.row(i).transpose());
    }
    return t_tcb;
  }

  Real LtToTcl(Real t_lt) {
    Real mjd_lt = TimeToMjd(t_lt);
    Real lt_tcl = -L_L / (1.0 - L_L) * (mjd_lt - MJD_COORDINATE_TT_TCG_TCB) * SECS_DAY;
    return t_lt - lt_tcl;
  }

  Real TclToLt(Real t_tcl) {
    Real mjd_tcl = TimeToMjd(t_tcl);
    Real lt_tcl = -L_L * (mjd_tcl - MJD_COORDINATE_TT_TCG_TCB) * SECS_DAY;
    return t_tcl + lt_tcl;
  }

  /**
   * @brief Calculate the proper time in Lunar Coordinate Time (TCL) frame
   * @param t_tcg Geocentric Coordinate Time (TCG)
   * @param x_mci Position vector in Moon-Centered Inertial (MCI) frame [m]
   * @return Proper time in TCL frame
   */
  Real GetProperTimeCorrectionTcl(Real t_tcg, const Vec3& x_mci) {
    Real t_tdb = TtToTdb(TcgToTt(t_tcg));
    Vec6 rv_LE = GetBodyPosVel(t_tdb, BodyId::EARTH, BodyId::MOON, Frame::GCRF);
    Vec3 r_vec_x_gcrf = ConvertFrame(t_tdb, x_mci, Frame::MOON_CI, Frame::GCRF);
    Vec3 r_vec_LE = rv_LE.head<3>();
    Vec3 v_vec_LE = rv_LE.tail<3>();
    Real proper_time_correction = -v_vec_LE.dot(r_vec_x_gcrf - r_vec_LE) / (C * C);
    return proper_time_correction;
  }

  VecX GetProperTimeCorrectionTcl(const VecX& t_tcg, const MatX& x_mci) {
    VecX proper_time_corrections(t_tcg.size());
    for (int i = 0; i < t_tcg.size(); i++) {
      proper_time_corrections(i) = GetProperTimeCorrectionTcl(t_tcg(i), x_mci.row(i).transpose());
    }
    return proper_time_corrections;
  }

  // =========================================================================
  // TL − TT conversion (Turyshev 2026, ApJ 997:97, Eq. 57)
  // =========================================================================
  //
  // Module-level optional Chebyshev fit of TL−TT as a function of TDB.
  // Populated by InitLtMinusTtFit(); cleared by ClearLtMinusTtFit().
  // When populated, TdbToLtMinusTt() uses it for fast evaluation; otherwise
  // it falls back to direct numerical integration from T_0.
  namespace {
    std::optional<ChebyshevFitModel> s_lt_minus_tt_fit;

  }

  /// @brief c⁻² integrand for TL−TT (Turyshev 2026, Eq. 57 inner integral).
  ///
  /// Returns  ½v²_EM + (GM_E − 2GM_M)/r_EM + W⊙_EM + L_G·(3/2)·GM_S/r_SE
  ///
  /// Integrated with sign −1/c² and added to the secular term gives TL−TT.
  static Real TdbToLtMinusTtIntegrandC2(Real t_tdb) {
    Vec6 rv_EM = GetBodyPosVel(t_tdb, BodyId::EARTH, BodyId::MOON, Frame::GCRF);
    Vec6 rv_ES = GetBodyPosVel(t_tdb, BodyId::SUN, BodyId::EARTH, Frame::GCRF);

    Vec3 r_EM = rv_EM.head<3>();
    Vec3 v_EM = rv_EM.tail<3>();
    Vec3 r_ES = rv_ES.head<3>();

    Real r_EM_norm = r_EM.norm();
    Real r_SE_norm = r_ES.norm();

    // v²_EM / 2
    Real kin = 0.5 * v_EM.squaredNorm();

    // (GM_E − 2GM_M) / r_EM  [note the factor of 2: −GM_M from TCL−TCG plus
    //  an additional −GM_M from the L_L scaling of TL = TCL − L_L·(TCL−T_0)]
    Real grav = (GM_EARTH - 2.0 * GM_MOON) / r_EM_norm;

    // Solar ℓ=2 tidal W⊙_EM = (3/2) GM_S / r_SE^5 · [(r_ES·r_EM)² − r_SE²r_EM²/3]
    Real dot_ESEM = r_ES.dot(r_EM);
    Real tidal = 1.5 * GM_SUN / pow(r_SE_norm, 5)
                 * (dot_ESEM * dot_ESEM - r_SE_norm * r_SE_norm * r_EM_norm * r_EM_norm / 3.0);

    // L_G correction: L_G · (3/2) · GM_S / r_SE
    Real sun_corr = L_G * 1.5 * GM_SUN / r_SE_norm;

    return kin + grav + tidal + sun_corr;
  }

  /// @brief c⁻⁴ integrand for TL−TT (Turyshev 2026, Eq. 57 second integral).
  ///
  /// Returns  3 · GM_S / r_SE · (v_E · v_EM)
  static Real TdbToLtMinusTtIntegrandC4(Real t_tdb) {
    Vec6 rv_EM = GetBodyPosVel(t_tdb, BodyId::EARTH, BodyId::MOON, Frame::GCRF);
    Vec6 rv_ES = GetBodyPosVel(t_tdb, BodyId::SUN, BodyId::EARTH, Frame::GCRF);
    Vec6 rv_E_bcrs = GetBodyPosVel(t_tdb, BodyId::SSB, BodyId::EARTH, Frame::GCRF);

    Vec3 v_EM = rv_EM.tail<3>();
    Vec3 v_E = rv_E_bcrs.tail<3>();
    Real r_SE_norm = rv_ES.head<3>().norm();

    return 3.0 * GM_SUN / r_SE_norm * v_E.dot(v_EM);
  }

  static Real TdbToLtMinusTtEarthMoonEndpoint(Real t_tdb) {
    Vec6 rv_EM = GetBodyPosVel(t_tdb, BodyId::EARTH, BodyId::MOON, Frame::GCRF);
    Vec6 rv_E_bcrs = GetBodyPosVel(t_tdb, BodyId::SSB, BodyId::EARTH, Frame::GCRF);
    return rv_E_bcrs.tail<3>().dot(rv_EM.head<3>());
  }

  /// @brief TL − TT as a function of TDB (Turyshev 2026, ApJ 997:97, Eq. 57).
  ///
  /// TL − TT = (L_G − L_L)/(1−L_B) · (TDB − T₀ − TDB₀)
  ///         − (1/c²) ∫_{T₀+TDB₀}^{TDB} [½v²_EM + (GM_E−2GM_M)/r_EM + W⊙_EM
  ///                                       + L_G·(3/2)·GM_S/r_SE] dTDB
  ///         + (1/c²)[v_E(T₀+TDB₀)·r_EM(T₀+TDB₀) − v_E(TDB)·r_EM(TDB)]
  ///         − (1/c⁴) ∫_{T₀+TDB₀}^{TDB} [3·GM_S/r_SE·(v_E·v_EM)] dTDB
  ///         + O(c⁻⁵)
  ///
  /// T₀  = coordinate time epoch (1977-01-01 00:00:32.184 TAI)
  ///     = MJD_COORDINATE_TT_TCG_TCB converted to seconds from J2000
  /// TDB₀ = −65.5 μs (DE405 TDB offset)
  /// T_L0 = T₀  (assumed → L_L·(T_L0−T₀) = 0)
  ///
  /// Integration uses a trapezoidal rule with dt = 0.01 days (≈ 864 s).
  /// If a Chebyshev fit has been initialised via InitLtMinusTtFit() and
  /// covers t_tdb, the fast fitted value is returned instead.
  ///
  /// @param t_tdb  Epoch in TDB [s from J2000]
  /// @return TL − TT [s]
  Real TdbToLtMinusTt(Real t_tdb) {
    // Fast path: use pre-fitted Chebyshev model if available
    {
      VecX val(1);
      bool ok = false;
#pragma omp critical
      {
        ok = s_lt_minus_tt_fit.has_value() && s_lt_minus_tt_fit->Eval(t_tdb, &val, nullptr);
      }
      if (ok) return val(0);
    }

    // Slow path: numerical integration from T₀ + TDB₀
    const double t0_s
        = (MJD_COORDINATE_TT_TCG_TCB - MJD_J2000_TT) * SECS_DAY + TDB_0;  // [s from J2000]
    const double t_end = t_tdb.val();

    if (t_end <= t0_s) {
      Real endpoint
          = (TdbToLtMinusTtEarthMoonEndpoint(t0_s) - TdbToLtMinusTtEarthMoonEndpoint(t_tdb))
            / (C * C);
      return (L_G - L_L) / (1.0 - L_B) * (t_tdb - t0_s) + endpoint;
    }

    const double dt_nom = 0.01 * SECS_DAY;  // 864 s nominal step
    double t_cur = t0_s;
    double int_c2 = 0.0, int_c4 = 0.0;
    double f_c2_prev = TdbToLtMinusTtIntegrandC2(t_cur).val();
    double f_c4_prev = TdbToLtMinusTtIntegrandC4(t_cur).val();

    while (t_cur < t_end) {
      double t_next = std::min(t_cur + dt_nom, t_end);
      double dt = t_next - t_cur;

      double f_c2 = TdbToLtMinusTtIntegrandC2(t_next).val();
      double f_c4 = TdbToLtMinusTtIntegrandC4(t_next).val();

      int_c2 += 0.5 * (f_c2 + f_c2_prev) * dt;
      int_c4 += 0.5 * (f_c4 + f_c4_prev) * dt;

      f_c2_prev = f_c2;
      f_c4_prev = f_c4;
      t_cur = t_next;
    }

    const double C2 = C * C;
    const double C4 = C2 * C2;
    Real secular = (L_G - L_L) / (1.0 - L_B) * (t_tdb - t0_s);
    Real endpoint
        = (TdbToLtMinusTtEarthMoonEndpoint(t0_s) - TdbToLtMinusTtEarthMoonEndpoint(t_tdb)) / C2;
    return secular - int_c2 / C2 - int_c4 / C4 + endpoint;
  }

  /// @brief Convert TDB epoch to Lunar Time (TL).
  /// @param t_tdb  Epoch in TDB [s from J2000]
  /// @return Epoch in TL [s from J2000]
  Real TdbToLt(Real t_tdb) { return TDBToTt(t_tdb) + TdbToLtMinusTt(t_tdb); }

  /// @brief Convert Lunar Time (TL) to TDB via Newton-Raphson inversion of TdbToLt.
  /// @param t_lt  Epoch in TL [s from J2000]
  /// @return Epoch in TDB [s from J2000]
  Real LtToTdb(Real t_lt) {
    Real t_tdb = t_lt;  // initial guess
    for (int iter = 0; iter < 10; iter++) {
      Real t_lt_est = TdbToLt(t_tdb);
      Real delta = t_lt - t_lt_est;
      t_tdb = t_tdb + delta;
      if (std::abs(delta.val()) < 1e-12) break;
    }
    return t_tdb;
  }

  /// @brief Convert TT epoch to Lunar Time (TL) via TDB.
  /// @param t_tt  Epoch in TT [s from J2000]
  /// @return Epoch in TL [s from J2000]
  Real TtToLt(Real t_tt) { return TdbToLt(TtToTdb(t_tt)); }

  /// @brief Convert Lunar Time (TL) to TT via inversion.
  /// @param t_lt  Epoch in TL [s from J2000]
  /// @return Epoch in TT [s from J2000]
  Real LtToTt(Real t_lt) { return TDBToTt(LtToTdb(t_lt)); }

  /// @brief Fit TL−TT(TDB) over a time window using piecewise Chebyshev polynomials.
  ///
  /// Performs a single forward sweep from the reference epoch T₀+TDB₀ to
  /// t_end_tdb, evaluating the Turyshev 2026 Eq. 57 integrand at every
  /// Chebyshev-Gauss node (processed in temporal order) and building the
  /// DCT Chebyshev fit segment-by-segment.  This avoids repeating the
  /// long T₀→t_start integral for each node independently.
  ///
  /// After this call, TdbToLtMinusTt() will use the fit for epochs in
  /// [t_start_tdb, t_end_tdb]; epochs outside that window fall back to
  /// direct integration.
  ///
  /// @param t_start_tdb    Start of fitting window [s from J2000, TDB]
  /// @param t_end_tdb      End   of fitting window [s from J2000, TDB]
  /// @param segment_length Segment length [s] (default 1 day)
  /// @param num_coeffs     Chebyshev degree+1 per segment (default 13)
  /// @brief Fit TDB-TCL over [t_start_tdb, t_end_tdb] with piecewise Chebyshev
  /// polynomials.
  ///
  /// The underlying TdbMinusTcl() integrates from T_0 (1977) on every call, so
  /// evaluating it once per Chebyshev node would cost O(nodes x span). All node
  /// times are collected, sorted, and covered by a SINGLE forward sweep instead,
  /// exactly as InitLtMinusTtFit does.
  void InitTdbMinusTclFit(Real t_start_tdb, Real t_end_tdb, double segment_length, int num_coeffs) {
    const double t_start = t_start_tdb.val();
    const double t_end = t_end_tdb.val();
    if (!(t_end > t_start)) throw std::runtime_error("InitTdbMinusTclFit: t_end must be > t_start");

    const Real t0_r = MjdToTime(MJD_COORDINATE_TT_TCG_TCB);
    const double t0_tcb = t0_r.val();

    // The integral is defined to be zero at T_0, so the sweep must START there.
    // Beginning earlier and zeroing the accumulator at that point would fold the
    // integral over [sweep_start, T_0] into every fitted value as a constant
    // offset (~3.7 s for a window reaching back to 1969).
    LUPNT_CHECK(t_start >= t0_tcb - 1.0,
                "InitTdbMinusTclFit: window must not start before T_0 (1977)",
                "InitTdbMinusTclFit");

    const double span = t_end - t_start;
    const int num_segs = std::max(1, static_cast<int>(std::ceil(span / segment_length)));
    const double seg_len = span / num_segs;

    struct NodeSpec {
      double t;  // TDB
      int seg_idx;
      int k;
    };
    std::vector<NodeSpec> nodes;
    nodes.reserve(static_cast<size_t>(num_segs) * num_coeffs);
    for (int sg = 0; sg < num_segs; sg++) {
      const double a = t_start + sg * seg_len;
      const double b = (sg == num_segs - 1) ? t_end : a + seg_len;
      const double mid = 0.5 * (a + b), rad = 0.5 * (b - a);
      for (int k = 0; k < num_coeffs; k++)
        nodes.push_back({mid + rad * std::cos((2.0 * k + 1.0) * PI / (2.0 * num_coeffs)), sg, k});
    }
    std::sort(nodes.begin(), nodes.end(),
              [](const NodeSpec& x, const NodeSpec& y) { return x.t < y.t; });

    // Single forward sweep of the TCB->TCL integrand, in TCB.
    const double dt_nom = 0.01 * SECS_DAY;
    const double C2 = C * C, C4 = C2 * C2;
    std::vector<std::vector<double>> vals(num_segs, std::vector<double>(num_coeffs, 0.0));

    auto tdb_to_tcb = [&](double t_tdb) { return (t_tdb - L_B * t0_tcb - TDB_0) / (1.0 - L_B); };

    double t_cur = t0_tcb;  // integral is zero here, by definition
    double i2 = 0.0, i4 = 0.0;
    double f2_prev = TcbToTclIntegrandC2(Real(t_cur)).val();
    double f4_prev = TcbToTclIntegrandC4(Real(t_cur)).val();

    size_t idx = 0;
    while (idx < nodes.size()) {
      const double t_node_tcb = tdb_to_tcb(nodes[idx].t);
      while (t_cur < t_node_tcb) {
        const double t_step = std::min(t_cur + dt_nom, t_node_tcb);
        const double dt = t_step - t_cur;
        const double f2 = TcbToTclIntegrandC2(Real(t_step)).val();
        const double f4 = TcbToTclIntegrandC4(Real(t_step)).val();
        i2 += 0.5 * (f2 + f2_prev) * dt;
        i4 += 0.5 * (f4 + f4_prev) * dt;
        f2_prev = f2;
        f4_prev = f4;
        t_cur = t_step;
      }
      const double val = -L_B * (t_cur - t0_tcb) + TDB_0 + i2 / C2 + i4 / C4;
      while (idx < nodes.size() && tdb_to_tcb(nodes[idx].t) <= t_cur) {
        vals[nodes[idx].seg_idx][nodes[idx].k] = val;
        idx++;
      }
    }

    ChebyshevFitModel model;
    model.t_start = t_start;
    model.t_end = t_end;
    model.num_dims = 1;
    model.num_coeffs = num_coeffs;
    model.segments.reserve(num_segs);
    for (int sg = 0; sg < num_segs; sg++) {
      const double a = t_start + sg * seg_len;
      const double b = (sg == num_segs - 1) ? t_end : a + seg_len;
      ChebyshevFitSegment seg;
      seg.t_mid = 0.5 * (a + b);
      seg.t_radius = 0.5 * (b - a);
      seg.coeffs.assign(1, std::vector<double>(num_coeffs, 0.0));
      for (int j = 0; j < num_coeffs; j++) {
        const double w = (j == 0) ? 1.0 : 2.0;
        for (int k = 0; k < num_coeffs; k++)
          seg.coeffs[0][j] += vals[sg][k] * std::cos((2.0 * k + 1.0) * PI * j / (2.0 * num_coeffs));
        seg.coeffs[0][j] *= w / num_coeffs;
      }
      model.segments.push_back(std::move(seg));
    }
#pragma omp critical
    {
      s_tdb_minus_tcl_fit = std::move(model);
    }
  }

  void ClearTdbMinusTclFit() {
#pragma omp critical
    {
      s_tdb_minus_tcl_fit.reset();
    }
  }

  bool HasFittedTdbMinusTcl(Real t_tdb) {
    bool ok = false;
#pragma omp critical
    {
      ok = s_tdb_minus_tcl_fit.has_value() && t_tdb.val() >= s_tdb_minus_tcl_fit->t_start
           && t_tdb.val() <= s_tdb_minus_tcl_fit->t_end;
    }
    return ok;
  }

  void SetTdbTclAutoFit(bool enable) {
#pragma omp critical
    {
      s_tdb_tcl_autofit = enable;
      if (enable) s_tdb_tcl_autofit_unavailable = false;
    }
  }

  bool GetTdbTclAutoFit() {
    bool v;
#pragma omp critical
    {
      v = s_tdb_tcl_autofit;
    }
    return v;
  }

  void InitLtMinusTtFit(Real t_start_tdb, Real t_end_tdb, double segment_length, int num_coeffs) {
    const double t_start = t_start_tdb.val();
    const double t_end = t_end_tdb.val();

    if (!(t_end > t_start)) throw std::runtime_error("InitLtMinusTtFit: t_end must be > t_start");

    const double t0_s = (MJD_COORDINATE_TT_TCG_TCB - MJD_J2000_TT) * SECS_DAY + TDB_0;

    // The LT-TT integral is defined to be zero at T_0, so the sweep must START
    // there. Beginning earlier (as std::min(t0_s, t_start) once did) and zeroing
    // the accumulator at that point folds the integral over [sweep_start, T_0]
    // into every fitted value as a constant offset -- the same defect that
    // InitTdbMinusTclFit guards against.
    LUPNT_CHECK(t_start >= t0_s - 1.0, "InitLtMinusTtFit: window must not start before T_0 (1977)",
                "InitLtMinusTtFit");
    const double t_sweep_start = t0_s;

    double span = t_end - t_start;
    int num_segs = std::max(1, static_cast<int>(std::ceil(span / segment_length)));
    double seg_len = span / num_segs;

    // -----------------------------------------------------------------------
    // 1. Collect ALL Chebyshev-Gauss node times across all segments (sorted).
    // -----------------------------------------------------------------------
    struct NodeSpec {
      double t;
      int seg_idx;
      int k;  // DCT index
    };
    std::vector<NodeSpec> all_nodes;
    all_nodes.reserve(num_segs * num_coeffs);

    for (int s = 0; s < num_segs; s++) {
      double seg_start = t_start + s * seg_len;
      double seg_end = (s == num_segs - 1) ? t_end : seg_start + seg_len;
      double mid = 0.5 * (seg_start + seg_end);
      double radius = 0.5 * (seg_end - seg_start);
      for (int k = 0; k < num_coeffs; k++) {
        double theta = (2.0 * k + 1.0) * PI / (2.0 * num_coeffs);
        double t = mid + radius * std::cos(theta);
        all_nodes.push_back({t, s, k});
      }
    }
    std::sort(all_nodes.begin(), all_nodes.end(),
              [](const NodeSpec& a, const NodeSpec& b) { return a.t < b.t; });

    // -----------------------------------------------------------------------
    // 2. Single forward sweep from t_sweep_start to t_end.
    //    Record TL−TT at each node (node_vals[seg][k]).
    // -----------------------------------------------------------------------
    const double dt_nom = 0.01 * SECS_DAY;  // 864 s integration step

    std::vector<std::vector<double>> node_vals(num_segs, std::vector<double>(num_coeffs, 0.0));

    double t_cur = t_sweep_start;
    double int_c2 = 0.0, int_c4 = 0.0;
    double f_c2_prev = TdbToLtMinusTtIntegrandC2(t_cur).val();
    double f_c4_prev = TdbToLtMinusTtIntegrandC4(t_cur).val();

    const double C2 = C * C;
    const double C4 = C2 * C2;

    auto lt_minus_tt_at = [&](double t_node) -> double {
      double secular = (L_G - L_L) / (1.0 - L_B) * (t_node - t0_s);
      double endpoint = (TdbToLtMinusTtEarthMoonEndpoint(t0_s).val()
                         - TdbToLtMinusTtEarthMoonEndpoint(t_node).val())
                        / C2;
      return secular - int_c2 / C2 - int_c4 / C4 + endpoint;
    };

    int node_idx = 0;
    int total_nodes = static_cast<int>(all_nodes.size());

    while (node_idx < total_nodes) {
      double t_node = all_nodes[node_idx].t;

      // Advance integration to t_node
      while (t_cur < t_node) {
        double t_step = std::min(t_cur + dt_nom, t_node);
        double dt = t_step - t_cur;

        double f_c2 = TdbToLtMinusTtIntegrandC2(t_step).val();
        double f_c4 = TdbToLtMinusTtIntegrandC4(t_step).val();

        int_c2 += 0.5 * (f_c2 + f_c2_prev) * dt;
        int_c4 += 0.5 * (f_c4 + f_c4_prev) * dt;

        f_c2_prev = f_c2;
        f_c4_prev = f_c4;
        t_cur = t_step;
      }

      // Record all nodes at exactly t_cur
      double val = lt_minus_tt_at(t_cur);
      while (node_idx < total_nodes && all_nodes[node_idx].t <= t_cur) {
        node_vals[all_nodes[node_idx].seg_idx][all_nodes[node_idx].k] = val;
        node_idx++;
      }
    }

    // -----------------------------------------------------------------------
    // 3. Build Chebyshev fit model from node values via DCT.
    // -----------------------------------------------------------------------
    ChebyshevFitModel model;
    model.t_start = t_start;
    model.t_end = t_end;
    model.num_dims = 1;
    model.num_coeffs = num_coeffs;
    model.segments.reserve(num_segs);

    for (int s = 0; s < num_segs; s++) {
      double seg_start = t_start + s * seg_len;
      double seg_end = (s == num_segs - 1) ? t_end : seg_start + seg_len;
      double mid = 0.5 * (seg_start + seg_end);
      double radius = 0.5 * (seg_end - seg_start);

      ChebyshevFitSegment seg;
      seg.t_mid = mid;
      seg.t_radius = radius;
      seg.coeffs.assign(1, std::vector<double>(num_coeffs, 0.0));

      for (int j = 0; j < num_coeffs; j++) {
        double w = (j == 0) ? 1.0 : 2.0;
        for (int k = 0; k < num_coeffs; k++) {
          double c = std::cos((2.0 * k + 1.0) * PI * j / (2.0 * num_coeffs));
          seg.coeffs[0][j] += node_vals[s][k] * c;
        }
        seg.coeffs[0][j] *= w / num_coeffs;
      }
      model.segments.push_back(std::move(seg));
    }

#pragma omp critical
    {
      s_lt_minus_tt_fit = std::move(model);
    }
  }

  /// @brief Clear the TL−TT Chebyshev fit, reverting to direct integration.
  void ClearLtMinusTtFit() {
#pragma omp critical
    {
      s_lt_minus_tt_fit.reset();
    }
  }

  /// @brief Return true if the TL−TT Chebyshev fit covers epoch t_tdb.
  bool HasFittedLtMinusTt(Real t_tdb) {
    bool ok = false;
#pragma omp critical
    {
      ok = s_lt_minus_tt_fit.has_value() && t_tdb.val() >= s_lt_minus_tt_fit->t_start
           && t_tdb.val() <= s_lt_minus_tt_fit->t_end;
    }
    return ok;
  }

  // =========================================================================
  // TT − TDB Chebyshev fit (high-fidelity DE440t path)
  // =========================================================================

  /// @brief Fit TT−TDB(TDB) over [t_start_tdb, t_end_tdb] with piecewise
  /// Chebyshev polynomials, sampling the high-fidelity DE440t TT-TDB ephemeris
  /// via spice::ConvertTime.
  ///
  /// After this call, TDBToTt() and TtToTdb() evaluate the fit directly (no
  /// SPICE call) for epochs inside the window, matching the DE440t ephemeris
  /// to sub-nanosecond accuracy while remaining safe to call from OpenMP
  /// parallel regions. Epochs outside the window fall back to the analytic
  /// series. This replaces per-call spice::ConvertTime, whose omp-critical
  /// SPICE access would otherwise serialize parallel time conversions.
  ///
  /// @param t_start_tdb    Start of the fitting window [s from J2000, TDB]
  /// @param t_end_tdb      End   of the fitting window [s from J2000, TDB]
  /// @param segment_length Segment length [s] (default 16 days)
  /// @param num_coeffs     Chebyshev degree+1 per segment (default 13)
  void InitTtMinusTdbFit(Real t_start_tdb, Real t_end_tdb, double segment_length, int num_coeffs) {
    auto sample = [](double t_tdb) -> VecXd {
      VecXd v(1);
      // Read the DE440t TT-TDB segment as a *clean offset*. Sampling via
      // spice::ConvertTime would return an absolute epoch, and subtracting
      // t_tdb back off would quantise every training node to the epoch ULP
      // (~245 ns at 2030) -- capping the fitted model's accuracy at that
      // level no matter how many Chebyshev coefficients are used.
      v(0) = spice::GetTimeEphemerisOffset(t_tdb, spice::kNaifTtMinusTdb).val();
      return v;
    };
    ChebyshevFitModel model = FitChebyshevModel(sample, t_start_tdb.val(), t_end_tdb.val(),
                                                /*num_dims=*/1, segment_length, num_coeffs);
#pragma omp critical
    {
      s_tt_minus_tdb_fit = std::move(model);
    }
  }

  /// @brief Clear the TT−TDB Chebyshev fit, reverting to the analytic series.
  void ClearTtMinusTdbFit() {
#pragma omp critical
    {
      s_tt_minus_tdb_fit.reset();
    }
  }

  /// @brief Return true if the TT−TDB Chebyshev fit covers epoch t.
  bool HasFittedTtMinusTdb(Real t) {
    bool ok = false;
#pragma omp critical
    {
      ok = s_tt_minus_tdb_fit.has_value() && t.val() >= s_tt_minus_tdb_fit->t_start
           && t.val() <= s_tt_minus_tdb_fit->t_end;
    }
    return ok;
  }

  VEC_IMP_REAL(UtcToUt1)
  VEC_IMP_REAL(Ut1ToUtc)
  VEC_IMP_REAL(TaiToUtc)
  VEC_IMP_REAL(UtcToTai)
  VEC_IMP_REAL(TaiToTt)
  VEC_IMP_REAL(TtToTai)
  VEC_IMP_REAL(TcgToTt)
  VEC_IMP_REAL(TtToTcg)
  VEC_IMP_REAL(TtToTdb)
  VEC_IMP_REAL(TDBToTt)
  VEC_IMP_REAL(TtMinusTdb)
  VEC_IMP_REAL(TdbMinusTcl)
  VEC_IMP_REAL(TdbMinusLt)
  VEC_IMP_REAL(TaiToGps)
  VEC_IMP_REAL(GpsToTai)
  VEC_IMP_REAL(TcbToTdb)
  VEC_IMP_REAL(TtToTcb)

  VEC_IMP_REAL(MjdToTime)
  VEC_IMP_REAL(TimeToMjd)
  VEC_IMP_REAL(JdToTime)
  VEC_IMP_REAL(TimeToJd)
  VEC_IMP_REAL(TclToLt)
  VEC_IMP_REAL(LtToTcl)
  VEC_IMP_REAL(TdbToLtMinusTt)
  VEC_IMP_REAL(TdbToLt)
  VEC_IMP_REAL(LtToTdb)
  VEC_IMP_REAL(TtToLt)
  VEC_IMP_REAL(LtToTt)

}  // namespace lupnt
