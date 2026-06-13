/**
 *  @file time_utils.cpp
 *  @author Keidai Iiyama
 *  @brief This file contains the implementation of time-related utility
 * functions.
 *  @version 0.1
 *  @date 2025-02-17
 */

#include "lupnt/environment/plasma/env/time_utils.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "lupnt/environment/plasma/core/math_utils.h"

namespace pecsim {

  void doy_to_mmdd(int year, int doy, int& mm, int& dd) {
    if (doy < 1 || doy > 366) {
      throw std::runtime_error("Day of year must be between 1 and 366.");
    }

    int month_days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    // Proper leap year check
    bool is_leap = (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
    if (is_leap) {
      month_days[1] = 29;
    }

    int i = 0;
    while (i < 12 && doy > month_days[i]) {
      doy -= month_days[i];
      ++i;
    }

    if (i >= 12) {
      throw std::runtime_error("Day of year exceeds total days in year.");
    }

    mm = i + 1;
    dd = doy;
  }

  void mmdd_to_doy(int year, int mm, int dd, int& doy) {
    if (mm < 1 || mm > 12) {
      throw std::runtime_error("Invalid month.");
    }

    int month_days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    bool is_leap = (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
    if (is_leap) {
      month_days[1] = 29;
    }

    if (dd < 1 || dd > month_days[mm - 1]) {
      throw std::runtime_error("Invalid day for the given month.");
    }

    doy = dd;
    for (int i = 0; i < mm - 1; ++i) {
      doy += month_days[i];
    }
  }

  std::array<int, 2> datetime_to_itime(DateTime datetime) {
    std::array<int, 2> itime;
    itime[0] = datetime.year * 1000 + datetime.doy;  // Year and day of the year
    int ihr = datetime.hour, imin = datetime.min, isec = datetime.sec;
    itime[1] = (ihr * 3600 + imin * 60 + isec) * 1000;  // Milliseconds of the day
    return itime;
  }

  DateTime itime_to_datetime(const std::array<int, 2>& itime) {
    DateTime datetime;
    datetime.year = itime[0] / 1000;  // Year
    datetime.doy = itime[0] % 1000;   // Day of the year

    int total_seconds = itime[1] / 1000;  // Convert milliseconds to seconds
    datetime.hour = total_seconds / 3600;
    total_seconds %= 3600;
    datetime.min = total_seconds / 60;
    datetime.sec = total_seconds % 60;

    return datetime;
  }

  double datetime_to_mjd(const DateTime& dt) {
    // Convert (year + doy) to month and day
    static const int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int month = 1, day = dt.doy;
    bool leap = (dt.year % 4 == 0 && (dt.year % 100 != 0 || dt.year % 400 == 0));
    int dim[12];
    std::copy(std::begin(days_in_month), std::end(days_in_month), dim);
    if (leap) dim[1] = 29;

    for (int i = 0; i < 12 && day > dim[i]; ++i) {
      day -= dim[i];
      ++month;
    }

    // Julian Date calculation
    int Y = dt.year;
    int M = month;
    if (M <= 2) {
      Y -= 1;
      M += 12;
    }
    int A = Y / 100;
    int B = 2 - A + A / 4;

    double D = day + (dt.hour + (dt.min + dt.sec / 60.0) / 60.0) / 24.0;

    double JD = std::floor(365.25 * (Y + 4716)) + std::floor(30.6001 * (M + 1)) + D + B - 1524.5;

    return JD - 2400000.5;
  }

  DateTime mjd_to_datetime(double mjd) {
    double jd = mjd + 2400000.5;
    double z = std::floor(jd + 0.5);
    double f = jd + 0.5 - z;

    double A = z;
    if (z >= 2299161) {
      int alpha = static_cast<int>((z - 1867216.25) / 36524.25);
      A = A + 1 + alpha - alpha / 4;
    }

    double B = A + 1524;
    int C = static_cast<int>((B - 122.1) / 365.25);
    int D = static_cast<int>(365.25 * C);
    int E = static_cast<int>((B - D) / 30.6001);

    double day = B - D - static_cast<int>(30.6001 * E) + f;
    int day_int = static_cast<int>(day);

    int month = (E < 14) ? E - 1 : E - 13;
    int year = (month > 2) ? C - 4716 : C - 4715;

    // Compute doy
    static const int days_before_month[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    int doy = days_before_month[month - 1] + day_int;
    if (month > 2 && (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0))) {
      doy += 1;  // leap year
    }

    double day_frac = day - day_int;
    int hour = static_cast<int>(day_frac * 24);
    int min = static_cast<int>((day_frac * 24 - hour) * 60);
    double sec = static_cast<double>(((day_frac * 24 - hour) * 60 - min) * 60);

    return DateTime{year, doy, hour, min, sec};
  }

  double long_to_lt(double along) {
    // Convert longitude in radians to magnetic local time (MLT) in hours
    double amlt = along / AMLTRAD;   // Convert radians to hours
    amlt = amod(amlt - 12.0, 24.0);  // Normalize to [0, 24) hours
    return amlt;
  }

  double lt_to_long(double amlt) {
    // Convert magnetic local time (MLT) to longitude in radians
    // MLT is in hours, so we convert it to radians
    double along = (amlt - 12.0) * AMLTRAD;
    return along;
  }

  double tj2000_to_mjd(double t_j2000) {
    // Convert TJD2000 to MJD
    double mjd = 51544.4996875 + t_j2000 / 86400.0;
    return mjd;
  }

  double mjd_to_tj2000(double mjd) {
    // Convert MJD to TJD2000
    double t_j2000 = (mjd - 51544.4996875) * 86400.0;
    return t_j2000;
  }

  double gregorian_to_mjd(int year, int month, int day, int hour, int min, double sec) {
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
    double frac_of_day = (hour + min / 60.0 + sec / 3600.0) / 24.0;
    return mjd_midnight + frac_of_day;
  }

}  // namespace pecsim
