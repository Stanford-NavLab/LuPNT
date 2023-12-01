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
  PA,        // Moon-Fixed with principal axes
  ME,        // Moon-Fixed with mean-Earth / polar axes
};

class CoordConverter {
 public:
  static Vector6 Convert(real epoch, Vector6 rv_in, CoordSystem coord_sys_in,
                         CoordSystem coord_sys_out);

 private:
  static Matrix6 ComputeITRFtoGCRF(real tai);
};
}  // namespace lupnt