/**
 * @file string_utils.h
 * @author Keidai Iiyama
 * @brief This file contains utility functions for string manipulation.
 * @version 0.1
 * @date 2025-02-17
 *
 */

#pragma once

namespace pecsim {

  const double RE = 6371.0;               // Earth radius in km
  const double PI = 3.1415927;            // Pi
  const double AMLTRAD = PI / 12.0;       // Radians per hour
  const double RAD2DEG = 180.0 / PI;      // Degrees per radian
  const double DEG2RAD = PI / 180.0;      // Radians per degree
  const double DELH = 1.0;                // Height step in RE units
  const double C = 299792.458;            // Speed of light [km/s]
  const double SECS_DAY = 86400.0;        // Seconds in a day
  const double TECU = 1e16;               // TEC unit in electrons/m^2
  const double GM_EARTH = 3.986004418e5;  // Gravitational constant for Earth [km^3/s^2]

}  // namespace pecsim
