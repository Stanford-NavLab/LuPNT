#include "lupnt/conversions/time_conversions.h"

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "lupnt/core/error.h"
#include "lupnt/environment/body.h"
#include "lupnt/environment/solar_system.h"
#include "lupnt/interfaces/eop.h"
#include "lupnt/interfaces/kernels.h"
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
      = {{Time::UT1, "UT1"},     {Time::UTC, "UTC"},      {Time::TAI, "TAI"}, {Time::TDB, "TDB"},
         {Time::TT, "TT"},       {Time::TCG, "TCG"},      {Time::TCB, "TCB"}, {Time::GPS, "GPS"},
         {Time::JD_TT, "JDTDT"}, {Time::JD_TDB, "JDTDB"}, {Time::TCL, "TCL"}, {Time::LT, "LT"}};

  const std::map<std::string, Time> string_to_time
      = {{"UT1", Time::UT1},     {"UTC", Time::UTC},      {"TAI", Time::TAI}, {"TDB", Time::TDB},
         {"TT", Time::TT},       {"TCG", Time::TCG},      {"TCB", Time::TCB}, {"GPS", Time::GPS},
         {"JDTDT", Time::JD_TT}, {"JDTDB", Time::JD_TDB}, {"TCL", Time::TCL}, {"LT", Time::LT}};

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
  Real ConvertTime(Real t, Time from, Time to) {
    if (from == to) return t;
    switch (from) {
      case Time::UT1: return ConvertTime(Ut1ToUtc(t), Time::UTC, to);
      case Time::UTC: {
        if (to == Time::UT1) return UtcToUt1(t);
        return ConvertTime(UtcToTai(t), Time::TAI, to);
      }
      case Time::TAI: {
        switch (to) {
          // If it reaches TAI, direct conversions to other systems are possible
          case Time::TAI: return t;
          case Time::GPS: return TaiToGps(t);
          case Time::UTC: return TaiToUtc(t);
          case Time::UT1: return UtcToUt1(TaiToUtc(t));
          case Time::TT: return TaiToTt(t);
          case Time::TCG: return TtToTcg(TaiToTt(t));
          case Time::TDB: return TtToTdb(TaiToTt(t));
          case Time::TCB: return TtToTcb(TaiToTt(t));
          case Time::TCL: return TcbToTcl(TtToTcb(TaiToTt(t)));
          default: break;
        }
        break;
      }
      case Time::TDB: {
        switch (to) {
          case Time::TT: return TDBToTt(t);
          case Time::TCB: return TtToTcb(TDBToTt(t));
          case Time::TCG: return TtToTcg(TDBToTt(t));
          case Time::TAI: return TtToTai(TDBToTt(t));
          case Time::LT: return TdbToLt(t);
          default: return ConvertTime(TtToTai(TDBToTt(t)), Time::TAI, to);
        }
      }
      case Time::TT: {
        switch (to) {
          case Time::TAI: return TtToTai(t);
          case Time::TDB: return TtToTdb(t);
          case Time::TCG: return TtToTcg(t);
          case Time::TCB: return TtToTcb(t);
          case Time::UTC: return TaiToUtc(TtToTai(t));
          case Time::UT1: return UtcToUt1(TaiToUtc(TtToTai(t)));
          case Time::GPS: return TaiToGps(TtToTai(t));
          case Time::LT: return TtToLt(t);
          default: return ConvertTime(TtToTai(t), Time::TAI, to);
        }
        break;
      }
      case Time::TCG: {
        switch (to) {
          default: return ConvertTime(TcgToTt(t), Time::TT, to);
        }
      }
      case Time::TCB: {
        if (to == Time::TCL) return TcbToTcl(t);
        return ConvertTime(TcbToTdb(t), Time::TDB, to);
      }
      case Time::GPS: return ConvertTime(GpsToTai(t), Time::TAI, to);
      case Time::TCL: {
        switch (to) {
          case Time::LT: return TclToLt(t);
          case Time::TCB: return TclToTcb(t);
          default: return ConvertTime(TclToTcb(t), Time::TCB, to);
        }
      }
      case Time::LT: {
        switch (to) {
          case Time::TCL: return LtToTcl(t);
          case Time::TDB: return LtToTdb(t);
          case Time::TT: return LtToTt(t);
          default: return ConvertTime(LtToTdb(t), Time::TDB, to);
        }
      }
    }
    throw std::invalid_argument("Invalid time conversion");
    return 0;
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

  /// @brief
  /// @param t_tdb
  /// @return
  /// @ref
  /// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/time.html#The%20Relationship%20between%20TT%20and%20TDB
  Real TDBToTt(Real t_tdb) {
    double k = 1.657e-3;
    double eb = 1.671e-2;
    Real mean_anom = 6.239996 + 1.99096871e-7 * t_tdb;
    Real ecc_anom = mean_anom + eb * sin(mean_anom);
    Real t_tt = t_tdb - k * sin(ecc_anom);
    return t_tt;
  }

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
  Real TtToTdb(Real t_tt) {
    Real days_j2000_tt = t_tt / SECS_DAY;
    Real me = (357.53 + 0.9856003 * days_j2000_tt) * RAD;
    Real t_tdb = t_tt + 0.001658 * sin(me) + 0.000014 * sin(2 * me);
    return t_tdb;
  }

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

    constexpr std::array<BodyId, 9> kTclExternalBodies = {
        BodyId::SUN,     BodyId::MERCURY, BodyId::VENUS,  BodyId::EARTH,   BodyId::MARS,
        BodyId::JUPITER, BodyId::SATURN,  BodyId::URANUS, BodyId::NEPTUNE,
    };

    Vec6 GetMoonBcrsStateForTcb(Real t_tcb) {
      return GetBodyPosVel(TcbToTdb(t_tcb), BodyId::SSB, BodyId::MOON, Frame::GCRF);
    }

    Real GetLunarExternalPotentialForTcb(Real t_tcb) {
      Real t_tdb = TcbToTdb(t_tcb);
      Vec3 x_moon = GetBodyPos(t_tdb, BodyId::SSB, BodyId::MOON, Frame::GCRF);
      Real w_ext = 0.0;

      for (BodyId body : kTclExternalBodies) {
        Vec3 x_body = GetBodyPos(t_tdb, BodyId::SSB, body, Frame::GCRF);
        w_ext += GetBodyGM(body) / (x_moon - x_body).norm();
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
      { ok = s_lt_minus_tt_fit.has_value() && s_lt_minus_tt_fit->Eval(t_tdb, &val, nullptr); }
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
  void InitLtMinusTtFit(Real t_start_tdb, Real t_end_tdb, double segment_length, int num_coeffs) {
    const double t_start = t_start_tdb.val();
    const double t_end = t_end_tdb.val();

    if (!(t_end > t_start)) throw std::runtime_error("InitLtMinusTtFit: t_end must be > t_start");

    const double t0_s = (MJD_COORDINATE_TT_TCG_TCB - MJD_J2000_TT) * SECS_DAY + TDB_0;
    const double t_sweep_start = std::min(t0_s, t_start);

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
    { s_lt_minus_tt_fit = std::move(model); }
  }

  /// @brief Clear the TL−TT Chebyshev fit, reverting to direct integration.
  void ClearLtMinusTtFit() {
#pragma omp critical
    { s_lt_minus_tt_fit.reset(); }
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
