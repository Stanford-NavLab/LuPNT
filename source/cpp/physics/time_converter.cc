#include "lupnt/physics/time_converter.h"

#include <functional>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <vector>

#include "lupnt/data/eop.h"
#include "lupnt/data/tai_utc.h"
#include "lupnt/numerics/graphs.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/solar_system.h"

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

  /// @brief Convert time from one time system to another
  /// @param t Time in the original time system
  /// @param from Original time system
  /// @param to Converted time system
  /// @return Real Time in the converted time system
  /// @note
  ///                     TCG
  ///                      |
  /// UT1 -- UTC -- TAI -- TT -> TCB
  ///                |     |      |
  ///               GPS   TDB <---+
  Real ConvertTime(Real t, Time from, Time to) {
    if (from == to) return t;
    switch (from) {
      case Time::UT1: return ConvertTime(UT12UTC(t), Time::UTC, to);
      case Time::UTC: {
        if (to == Time::UT1) return UTC2UT1(t);
        return ConvertTime(UTC2TAI(t), Time::TAI, to);
      }
      case Time::TAI: {
        switch (to) {
          case Time::TAI: return t;
          case Time::GPS: return TAI2GPS(t);
          case Time::UTC: return TAI2UTC(t);
          case Time::UT1: return UTC2UT1(TAI2UTC(t));
          case Time::TT: return TAI2TT(t);
          case Time::TCG: return TT2TCG(TAI2TT(t));
          case Time::TDB: return TT2TDB(TAI2TT(t));
          case Time::TCB: return TT2TCB(TAI2TT(t));
          case Time::JD_TT: return Time2JD(TAI2TT(t));
          case Time::JD_TDB: return Time2JD(TT2TDB(TAI2TT(t)));
          default: break;
        }
      }
      case Time::TDB: {
        switch (to) {
          case Time::TT: return TDB2TT(t);
          case Time::TCB: return TT2TCB(TDB2TT(t));
          case Time::TCG: return TT2TCG(TDB2TT(t));
          case Time::TAI: return TT2TAI(TDB2TT(t));
          case Time::JD_TDB: return Time2JD(t);
          case Time::JD_TT: return Time2JD(TDB2TT(t));
          default: return ConvertTime(TT2TAI(TDB2TT(t)), Time::TAI, to);
        }
      }
      case Time::TT: {
        switch (to) {
          case Time::TAI: return TT2TAI(t);
          case Time::TDB: return TT2TDB(t);
          case Time::TCG: return TT2TCG(t);
          case Time::TCB: return TT2TCB(t);
          case Time::UTC: return TAI2UTC(TT2TAI(t));
          case Time::UT1: return UTC2UT1(TAI2UTC(TT2TAI(t)));
          case Time::GPS: return TAI2GPS(TT2TAI(t));
          case Time::JD_TT: return Time2JD(t);
          case Time::JD_TDB: return Time2JD(TT2TDB(t));
          default: break;
        }
      }
      case Time::TCG: return ConvertTime(TCG2TT(t), Time::TT, to);
      case Time::TCB: return ConvertTime(TCB2TDB(t), Time::TDB, to);
      case Time::GPS: return ConvertTime(GPS2TAI(t), Time::TAI, to);
      default: break;
    }
    throw std::invalid_argument("Invalid time conversion");
    return 0;
  }

  VecX ConvertTime(VecX t, Time from, Time to) {
    VecX t_out(t.size());
    for (int i = 0; i < t.size(); i++) {
      t_out(i) = ConvertTime(t(i), from, to);
    }
    return t_out;
  }

  Real UTC2UT1(Real t_utc) {
    Real mjd_utc = t_utc / SECS_DAY + MJD_J2000;
    Real ut1_utc = GetUt1UtcDifference(mjd_utc);
    Real t_ut1 = t_utc + ut1_utc;
    return t_ut1;
  }

  Real UT12UTC(Real t_ut1) {
    Real mjd_ut1 = t_ut1 / SECS_DAY + MJD_J2000;
    Real ut1_utc = GetUt1UtcDifference(mjd_ut1);
    Real t_utc = t_ut1 - ut1_utc;
    return t_utc;
  }

  Real UTC2TAI(Real t_utc) {
    Real mjd_utc = t_utc / SECS_DAY + MJD_J2000;
    Real tai_utc = GetTaiUtcDifference(mjd_utc.val());
    Real t_tai = t_utc + tai_utc;
    return t_tai;
  }

  Real TAI2UTC(Real t_tai) {
    Real mjd_tai = t_tai / SECS_DAY + MJD_J2000;
    Real tai_utc = GetTaiUtcDifference(mjd_tai.val());
    Real t_utc = t_tai - tai_utc;
    return t_utc;
  }

  Real TAI2TT(Real t_tai) { return t_tai + TT_TAI_OFFSET; }

  Real TT2TAI(Real t_tt) { return t_tt - TT_TAI_OFFSET; }

  Real TT2TCG(Real t_tt) {
    Real jd_tt = JD_J2000 + t_tt / SECS_DAY;
    Real tt_tcg = -L_G / (1.0 - L_G) * (jd_tt - JD_T0) * SECS_DAY;
    return t_tt - tt_tcg;
  }

  Real TCG2TT(Real t_tcg) {
    Real jd_tcg = JD_J2000 + t_tcg / SECS_DAY;
    Real tt_tcg = -L_G * (jd_tcg - JD_T0) * SECS_DAY;
    return t_tcg + tt_tcg;
  }

  /// @brief
  /// @param t_tdb
  /// @return
  /// @ref
  /// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/time.html#The%20Relationship%20between%20TT%20and%20TDB
  Real TDB2TT(Real t_tdb) {
    double k = 1.657e-3;
    double eb = 1.671e-2;
    Real mean_anom = 6.239996 + 1.99096871e-7 * t_tdb;
    Real ecc_anom = mean_anom + eb * sin(mean_anom);
    Real t_tt = t_tdb - k * sin(ecc_anom);
    return t_tt;
  }

  /// @brief
  /// @param t_tt
  /// @return
  /// @ref
  /// https://gssc.esa.int/navipedia/index.php/Transformations_between_Time_Systems#TDT_-_TDB.2C_TCB
  /// @note Accurate to about 30 microseconds
  Real TT2TDB(Real t_tt) {
    double k = 1.657e-3;
    double eb = 1.671e-2;
    Real mean_anom = 6.239996 + 1.99096871e-7 * t_tt;
    Real ecc_anom = mean_anom + eb * sin(mean_anom);
    Real t_tdb = t_tt + k * sin(ecc_anom);
    return t_tdb;
  }

  Real TAI2GPS(Real t_tai) {
    Real t_gps = t_tai - 19.0;
    return t_gps;
  }

  Real GPS2TAI(Real t_gps) {
    Real t_tai = t_gps + 19.0;
    return t_tai;
  }

  Real TCB2TDB(Real t_tcb) {
    const double tdb0 = -6.44e-5;
    Real jd_tcb = JD_J2000 + t_tcb / SECS_DAY;
    Real t_tdb = 1.55051976772e-8 * (jd_tcb - JD_T0) * SECS_DAY + tdb0;
    return t_tdb;
  }

  Real TT2TCB(Real t_tt) {
    Real t_tai = TT2TAI(t_tt);
    Real jd_tai = JD_J2000 + t_tai / SECS_DAY;
    Real t_tcb = t_tt + 1.5505197677e-8 * (jd_tai - JD_T0) * SECS_DAY;
    return t_tcb;
  }

  Real EarthRotationAngle(Real t_ut1) {
    double theta_0 = 0.7790572732640;
    double dtheta_dt = 1.00273781191135448;
    Real theta_era = TWO_PI * (theta_0 + dtheta_dt * t_ut1 / SECS_DAY);
    return Wrap2Pi(theta_era);
  }

  Real Gregorian2MJD(int year, int month, int day, int hour, int min, Real sec) {
    if (month <= 2) {
      month += 12;
      --year;
    }
    int b;
    if ((10000L * year + 100L * month + day) <= 15821004L)
      b = -2 + ((year + 4716) / 4) - 1179;  // Julian calendar
    else
      b = (year / 400) - (year / 100) + (year / 4);  // Gregorian calendar

    Real mjd_midnight = 365L * year - 679004L + b + int(30.6001 * (month + 1)) + day;
    Real frac_of_day = (hour + min / 60.0 + sec / 3600.0) / 24.0;
    return mjd_midnight + frac_of_day;
  }

  std::tuple<int, int, int, int, int, Real> MJD2Gregorian(Real mjd) {
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

  Real Gregorian2Time(int year, int month, int day, int hour, int min, Real sec) {
    Real mjd = Gregorian2MJD(year, month, day, hour, min, sec);
    return MJD2Time(mjd);
  }

  /// @brief Greenwich Mean Sidereal Time
  /// @param mjd_ut1 UT1 (Modified Julian Date)
  /// @return GMST [rad]
  Real GreenwichMeanSiderealTime(Real mjd_ut1) {
    Real mjd0 = floor(mjd_ut1);
    Real ut1 = SECS_DAY * (mjd_ut1 - mjd0);  // [s]
    Real T0 = (mjd0 - MJD_J2000) / DAYS_CENTURY;
    Real T = (mjd_ut1 - MJD_J2000) / DAYS_CENTURY;

    Real gmst = 24110.54841 + 8640184.812866 * T0 + 1.002737909350795 * ut1
                + (0.093104 - 6.2e-6 * T) * T * T;  // [s]

    return TWO_PI * frac(gmst / SECS_DAY);  // [rad]
  }

  Real MJD2Time(Real mjd) { return (mjd - MJD_J2000) * SECS_DAY; }

  Real Time2MJD(Real t) { return t / SECS_DAY + MJD_J2000; }

  Real JD2Time(Real jd) { return (jd - JD_J2000) * SECS_DAY; }

  Real Time2JD(Real t) { return t / SECS_DAY + JD_J2000; }

  /// @brief Convert Modified Julian Date to date string
  /// @param mjd Modified Julian Date
  /// @param precision Number of seconds precision
  /// @return Date string
  std::string MJD2GregorianString(Real mjd, int precision) {
    double pow10 = pow(10, precision);
    Real mjd_round = (round(mjd * SECS_DAY * pow10, precision) + 0.1) / (SECS_DAY * pow10);
    auto [year, month, day, hour, min, sec] = MJD2Gregorian(mjd_round);
    std::stringstream ss;
    sec = round(sec, precision);
    ss << year << "/";
    ss << std::setw(2) << std::setfill('0') << month << "/";
    ss << std::setw(2) << std::setfill('0') << day << " ";
    ss << std::setw(2) << std::setfill('0') << hour << ":";
    ss << std::setw(2) << std::setfill('0') << min << ":";
    ss << std::setw(2) << std::setfill('0') << floor(sec) << ".";
    ss << std::fixed << std::setprecision(0) << std::setw(precision) << std::setfill('0')
       << round((sec - floor(sec)) * pow(10, precision));
    return ss.str();
  }

  std::string Time2GregorianString(Real t, int precision) {
    Real mjd = Time2MJD(t);
    return MJD2GregorianString(mjd, precision);
  }

  /// @brief Greenwich Apparent Sidereal Time
  /// @param mjd_ut1 UT1 (Modified Julian Date)
  /// @return GAST [rad]
  Real GreenwichApparentSiderealTime(Real mjd_ut1) {
    return mod(GreenwichMeanSiderealTime(mjd_ut1) + EquinoxEquation(mjd_ut1), TWO_PI);
  }

  VEC_IMP_REAL(UTC2UT1)
  VEC_IMP_REAL(UT12UTC)
  VEC_IMP_REAL(TAI2UTC)
  VEC_IMP_REAL(UTC2TAI)
  VEC_IMP_REAL(TAI2TT)
  VEC_IMP_REAL(TT2TAI)
  VEC_IMP_REAL(TCG2TT)
  VEC_IMP_REAL(TT2TCG)
  VEC_IMP_REAL(TT2TDB)
  VEC_IMP_REAL(TDB2TT)
  VEC_IMP_REAL(TAI2GPS)
  VEC_IMP_REAL(GPS2TAI)
  VEC_IMP_REAL(TCB2TDB)
  VEC_IMP_REAL(TT2TCB)

  VEC_IMP_REAL(MJD2Time)
  VEC_IMP_REAL(Time2MJD)
  VEC_IMP_REAL(JD2Time)
  VEC_IMP_REAL(Time2JD)

}  // namespace lupnt
