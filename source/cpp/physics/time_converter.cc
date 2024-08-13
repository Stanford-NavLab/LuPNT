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
  {{TimeSys::from, TimeSys::to}, [](Real t) -> Real { return func(t); }}

/// @ref
/// D. Folta, N. Bosanac, I. Elliott, L. Mann, R. Mesarch, and J. Rosales,
/// ‘Astrodynamics Convention and Modeling Reference for Lunar, Cislunar, and
/// Libration Point Orbits’, Jan. 2022.
/// @ref
// O. Montenbruck and G. Eberhard, “Satellite Orbits: Models, Methods, and
// Applications,” Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.

namespace lupnt {

  std::map<std::pair<std::string, std::string>, std::function<Real(Real)>> time_conversions
      = {TIME_CONVERSION(UTC, UT1, UTCtoUT1), TIME_CONVERSION(UT1, UTC, UT1toUTC),
         TIME_CONVERSION(TAI, UTC, TAItoUTC), TIME_CONVERSION(UTC, TAI, UTCtoTAI),
         TIME_CONVERSION(TAI, TT, TAItoTT),   TIME_CONVERSION(TT, TAI, TTtoTAI),
         TIME_CONVERSION(TCG, TT, TCGtoTT),   TIME_CONVERSION(TT, TCG, TTtoTCG),
         TIME_CONVERSION(TT, TDB, TTtoTDB),   TIME_CONVERSION(TDB, TT, TDBtoTT),
         TIME_CONVERSION(TAI, GPS, TAItoGPS), TIME_CONVERSION(GPS, TAI, GPStoTAI),
         TIME_CONVERSION(TCB, TDB, TCBtoTDB), TIME_CONVERSION(TT, TCB, TTtoTCB)};

  Real ConvertTime(Real t, const std::string& from, const std::string& to) {
    if (from == to) return t;
    std::vector<std::string> path = FindShortestPath(from, to, time_conversions);
    Real t_out = t;
    for (size_t i = 0; i < path.size() - 1; i++) {
      std::function<Real(Real)> f = time_conversions[{path[i], path[i + 1]}];
      t_out = f(t_out);
    }
    return t_out;
  }

  VecX ConvertTime(VecX t, const std::string& from, const std::string& to) {
    VecX t_out(t.size());
    for (int i = 0; i < t.size(); i++) {
      t_out(i) = ConvertTime(t(i), from, to);
    }
    return t_out;
  }

  Real UTCtoUT1(Real t_utc) {
    Real mjd_utc = t_utc / SECS_DAY + MJD_J2000;
    Real ut1_utc = GetUt1UtcDifference(mjd_utc);
    Real t_ut1 = t_utc + ut1_utc;
    return t_ut1;
  }

  Real UT1toUTC(Real t_ut1) {
    Real mjd_ut1 = t_ut1 / SECS_DAY + MJD_J2000;
    Real ut1_utc = GetUt1UtcDifference(mjd_ut1);
    Real t_utc = t_ut1 - ut1_utc;
    return t_utc;
  }

  Real UTCtoTAI(Real t_utc) {
    Real mjd_utc = t_utc / SECS_DAY + MJD_J2000;
    Real tai_utc = GetTaiUtcDifference(mjd_utc.val());
    Real t_tai = t_utc + tai_utc;
    return t_tai;
  }

  Real TAItoUTC(Real t_tai) {
    Real mjd_tai = t_tai / SECS_DAY + MJD_J2000;
    Real tai_utc = GetTaiUtcDifference(mjd_tai.val());
    Real t_utc = t_tai - tai_utc;
    return t_utc;
  }

  Real TAItoTT(Real t_tai) { return t_tai + TT_TAI_OFFSET; }

  Real TTtoTAI(Real t_tt) { return t_tt - TT_TAI_OFFSET; }

  Real TTtoTCG(Real t_tt) {
    Real jd_tt = JD_J2000 + t_tt / SECS_DAY;
    Real tt_tcg = -L_G / (1.0 - L_G) * (jd_tt - JD_T0) * SECS_DAY;
    return t_tt - tt_tcg;
  }

  Real TCGtoTT(Real t_tcg) {
    Real jd_tcg = JD_J2000 + t_tcg / SECS_DAY;
    Real tt_tcg = -L_G * (jd_tcg - JD_T0) * SECS_DAY;
    return t_tcg + tt_tcg;
  }

  /// @brief
  /// @param t_tdb
  /// @return
  /// @ref
  /// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/time.html#The%20Relationship%20between%20TT%20and%20TDB
  Real TDBtoTT(Real t_tdb) {
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
  Real TTtoTDB(Real t_tt) {
    double k = 1.657e-3;
    double eb = 1.671e-2;
    Real mean_anom = 6.239996 + 1.99096871e-7 * t_tt;
    Real ecc_anom = mean_anom + eb * sin(mean_anom);
    Real t_tdb = t_tt + k * sin(ecc_anom);
    return t_tdb;
  }

  Real TAItoGPS(Real t_tai) {
    Real t_gps = t_tai - 19.0;
    return t_gps;
  }

  Real GPStoTAI(Real t_gps) {
    Real t_tai = t_gps + 19.0;
    return t_tai;
  }

  Real TCBtoTDB(Real t_tcb) {
    const double tdb0 = -6.44e-5;
    Real jd_tcb = JD_J2000 + t_tcb / SECS_DAY;
    Real t_tdb = 1.55051976772e-8 * (jd_tcb - JD_T0) * SECS_DAY + tdb0;
    return t_tdb;
  }

  Real TTtoTCB(Real t_tt) {
    Real t_tai = TTtoTAI(t_tt);
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

  Real GregorianToMJD(int year, int month, int day, int hour, int min, Real sec) {
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

  std::tuple<int, int, int, int, int, Real> MJDtoGregorian(Real mjd) {
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
    Real mjd = GregorianToMJD(year, month, day, hour, min, sec);
    return MJDtoTime(mjd);
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

  Real MJDtoTime(Real mjd) { return (mjd - MJD_J2000) * SECS_DAY; }

  Real TimeToMJD(Real t) { return t / SECS_DAY + MJD_J2000; }

  Real JDtoTime(Real jd) { return (jd - JD_J2000) * SECS_DAY; }

  Real TimeToJD(Real t) { return t / SECS_DAY + JD_J2000; }

  /// @brief Convert Modified Julian Date to date string
  /// @param mjd Modified Julian Date
  /// @param precision Number of seconds precision
  /// @return Date string
  std::string MJDtoGregorianString(Real mjd, int precision) {
    double pow10 = pow(10, precision);
    Real mjd_round = (round(mjd * SECS_DAY * pow10, precision) + 0.1) / (SECS_DAY * pow10);
    auto [year, month, day, hour, min, sec] = MJDtoGregorian(mjd_round);
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

  std::string TimeToGregorianString(Real t, int precision) {
    Real mjd = TimeToMJD(t);
    return MJDtoGregorianString(mjd, precision);
  }

  /// @brief Greenwich Apparent Sidereal Time
  /// @param mjd_ut1 UT1 (Modified Julian Date)
  /// @return GAST [rad]
  Real GreenwichApparentSiderealTime(Real mjd_ut1) {
    return mod(GreenwichMeanSiderealTime(mjd_ut1) + EquinoxEquation(mjd_ut1), TWO_PI);
  }

  VEC_IMP_REAL(UTCtoUT1)
  VEC_IMP_REAL(UT1toUTC)
  VEC_IMP_REAL(TAItoUTC)
  VEC_IMP_REAL(UTCtoTAI)
  VEC_IMP_REAL(TAItoTT)
  VEC_IMP_REAL(TTtoTAI)
  VEC_IMP_REAL(TCGtoTT)
  VEC_IMP_REAL(TTtoTCG)
  VEC_IMP_REAL(TTtoTDB)
  VEC_IMP_REAL(TDBtoTT)
  VEC_IMP_REAL(TAItoGPS)
  VEC_IMP_REAL(GPStoTAI)
  VEC_IMP_REAL(TCBtoTDB)
  VEC_IMP_REAL(TTtoTCB)

  VEC_IMP_REAL(MJDtoTime)
  VEC_IMP_REAL(TimeToMJD)
  VEC_IMP_REAL(JDtoTime)
  VEC_IMP_REAL(TimeToJD)

}  // namespace lupnt
