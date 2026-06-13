/**
 *  @file env/time_utils.h
 *  @author Keidai Iiyama
 *  @brief This file contains the interface for time-related utility functions.
 *  @version 0.1
 *  @date 2025-02-17
 *
 */

#pragma once

#include <array>

namespace pecsim {

  struct DateTime {
    int year;
    int doy;
    int hour;
    int min;
    double sec;
  };

  /**
   * @brief Convert day of year to month and day.
   *
   * @param year The year.
   * @param doy The day of the year (1-366).
   * @param mm Output parameter for month (1-12).
   * @param dd Output parameter for day (1-31).
   */
  void doy_to_mmdd(int year, int doy, int& mm, int& dd);

  /**
   * @brief Convert month and day to day of year.
   *
   * @param year The year.
   * @param mm The month (1-12).
   * @param dd The day (1-31).
   * @param doy Output parameter for day of the year (1-366).
   */
  void mmdd_to_doy(int year, int mm, int dd, int& doy);

  /**
   * @brief Convert DateTime to itime array.
   *
   * @param datetime The DateTime object.
   * @return std::array<int, 2> The itime array.
   */
  std::array<int, 2> datetime_to_itime(DateTime datetime);

  /*
   * @brief Convert itime array to DateTime object.
   *
   * @param itime The itime array, where itime[0] is year and day of year,
   *              and itime[1] is milliseconds of the day.
   * @return DateTime The DateTime object.
   */
  DateTime itime_to_datetime(const std::array<int, 2>& itime);

  /**
   * @brief Convert DateTime to Modified Julian Date (MJD).
   *
   * @param datetime The DateTime object.
   * @return double The Modified Julian Date.
   */
  double datetime_to_mjd(const DateTime& datetime);

  /**
   * @brief Convert Modified Julian Date (MJD) to DateTime.
   *
   * @param mjd The Modified Julian Date.
   * @return DateTime The DateTime object.
   */
  DateTime mjd_to_datetime(double mjd);

  /**
   * @brief Convert longitude in degrees to local time in hours.
   *
   * @param along Longitude in radians.
   * @return double Local time in hours.
   */
  double long_to_lt(double along);

  /**
   * @brief Convert magnetic local time (MLT) in hours to longitude in radians.
   *
   * @param amlt Magnetic local time in hours.
   * @return double Longitude in radians.
   */
  double lt_to_long(double amlt);

  /**
   * @brief Convert TJD2000 to Modified Julian Date (MJD).
   *
   * @param t_j2000 Time (seconds since J2000 epoch).
   * @return double The Modified Julian Date.
   */
  double tj2000_to_mjd(double t_j2000);

  /**
   * @brief Convert Modified Julian Date (MJD) to TJD2000.
   *
   * @param mjd The Modified Julian Date.
   * @return double Time (seconds since J2000 epoch).
   */
  double mjd_to_tj2000(double mjd);

  /**
   * @brief Convert Gregorian date to Modified Julian Date (MJD).
   *
   * @param year The year.
   * @param month The month (1-12).
   * @param day The day (1-31).
   * @param hour The hour (0-23).
   * @param min The minute (0-59).
   * @param sec The second (0.0-59.999...).
   * @return double The Modified Julian Date.
   */
  double gregorian_to_mjd(int year, int month, int day, int hour, int min, double sec);

}  // namespace pecsim
