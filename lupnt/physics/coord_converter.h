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

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <memory>

#include "cheby.h"
#include "lupnt/core/constants.h"

namespace ad = autodiff;

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

  static CoordSystem GetCoordTypeID(const std::string &str);

  static VectorXreal Convert(const VectorXreal rv_in,
                                 const real epoch,
                                 const std::string coord_sys_in,
                                 const std::string coord_sys_out);
  static VectorXreal Convert(const VectorXreal rv_in,
                                 const real epoch,
                                 const CoordSystem coord_sys_in,
                                 const CoordSystem coord_sys_out);

 private:
  static MatrixXreal ComputeITRFtoGCRF(const real tai);
  static Matrix3real R1(real phi);
  static Matrix3real R2(real phi);
  static Matrix3real R3(real phi);
  static Matrix3real Skew(Vector3real vec);
};
}  // namespace lupnt