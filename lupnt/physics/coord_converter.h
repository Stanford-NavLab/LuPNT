/**
 * @file CoordConverter.cpp
 * @author Stanford NAV LAB
 * @brief Coordinate conversion functions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <memory>

#include "cheby.h"
#include "lupnt/core/constants.h"

namespace lupnt {

enum CoordSystem {
  ITRF = 0,  // International Terrestrial Reference Frame
  GCRF,      // Geocentric Reference System
  ICRF,      // International Celestial Reference System
  SER,       // Sun-Earth Rotating Frame
  GSE,       // Geocentric Solar Ecliptic
  EME,       // Earth-Centered mean equator and equinox at J2000 epoch
  MOD,       // Mean of date equatorial system
  TOD,       // True of date equatorial system
  EMR,       // Earth-Moon Rotating Frame
  MI,        // Moon-centered Inertial Frame  (Axis aligened with ICRF)
  PA,        // Moon Fixed with principal axes
  ME,        // Moon-Fixed with mean-Earth / polar axes
  RTN,       // Radial-Tangential-Normal
  CoordSystemCount,
  NONE,
};

class CoordConverter {
 public:
  static const std::string COORD_SYSTEM_TEXT[CoordSystemCount];

  static Vector6real Convert(Vector6real rv_in, real epoch,
                             CoordSystem coord_sys_in,
                             CoordSystem coord_sys_out);

 private:
  static Matrix6real ComputeITRFtoGCRF(real tai);
  static Matrix3real R1(real phi);
  static Matrix3real R2(real phi);
  static Matrix3real R3(real phi);
  static Matrix3real Skew(Vector3real x);
};
}  // namespace lupnt