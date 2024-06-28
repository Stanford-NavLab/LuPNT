#include "time.h"

#include "lupnt/numerics/math_utils.h"

namespace lupnt {

/// @ref
// O. Montenbruck and G. Eberhard, “Satellite Orbits: Models, Methods, and
// Applications,” Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.
Real Calendar2ModJulianDate(int year, int month, int day, int hour, int min,
                            Real sec) {
  if (month <= 2) {
    month += 12;
    --year;
  }

  int b;
  if ((10000L * year + 100L * month + day) <= 15821004L)
    b = -2 + ((year + 4716) / 4) - 1179;  // Julian calendar
  else
    b = (year / 400) - (year / 100) + (year / 4);  // Gregorian calendar

  Real mjd_midnight =
      365L * year - 679004L + b + int(30.6001 * (month + 1)) + day;
  Real frac_of_day = (hour + min / 60.0 + sec / 3600.0) / 24.0;

  return mjd_midnight + frac_of_day;
}

/// @ref
// O. Montenbruck and G. Eberhard, “Satellite Orbits: Models, Methods, and
// Applications,” Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.
std::tuple<int, int, int, int, int, Real> ModJulianDate2Calendar(Real mjd) {
  long a, b, c, d, e, f;

  // Convert Julian day number to calendar date
  a = long(mjd + 2400001.0);

  if (a < 2299161) {  // Julian calendar
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
  Real x = (hours - hour) * MINUTES_HOUR;
  int min = int(x);
  Real sec = (x - min) * SECS_MINUTE;

  return std::make_tuple(year, month, day, hour, min, sec);
}

/// @ref
// O. Montenbruck and G. Eberhard, “Satellite Orbits: Models, Methods, and
// Applications,” Berlin : New York: Springer, 2000.
// doi: 10.1007/978-3-642-58351-3.
Real GreenwichMeanSiderealTime(Real mjd_ut1) {
  Real mjd0 = floor(mjd_ut1);
  Real ut1 = SECS_DAY * (mjd_ut1 - mjd0);  // [s]
  Real T0 = (mjd0 - MJD_J2000) / DAYS_JULIAN_CENTURY;
  Real T = (mjd_ut1 - MJD_J2000) / DAYS_JULIAN_CENTURY;

  Real gmst = 24110.54841 + 8640184.812866 * T0 + 1.002737909350795 * ut1 +
              (0.093104 - 6.2e-6 * T) * T * T;  // [s]

  return TWO_PI * frac(gmst / SECS_DAY);  // [rad]
}

/// @brief Convert Modified Julian Date to date string
/// @param mjd Modified Julian Date
/// @param precision Number of seconds precision
/// @return Date string
std::string FormatDate(Real mjd, int precision) {
  double pow10 = pow(10, precision);
  Real mjd_round =
      (round(mjd * SECS_DAY * pow10, precision) + 0.1) / (SECS_DAY * pow10);
  auto [year, month, day, hour, min, sec] = ModJulianDate2Calendar(mjd_round);

  std::stringstream ss;
  sec = round(sec, precision);
  ss << year << "/";
  ss << std::setw(2) << std::setfill('0') << month << "/";
  ss << std::setw(2) << std::setfill('0') << day << " ";
  ss << std::setw(2) << std::setfill('0') << hour << ":";
  ss << std::setw(2) << std::setfill('0') << min << ":";
  ss << std::setw(2) << std::setfill('0') << floor(sec) << ".";
  ss << std::setw(precision) << std::setfill('0')
     << round((sec - floor(sec)) * pow(10, precision));
  return ss.str();
}
}  // namespace lupnt