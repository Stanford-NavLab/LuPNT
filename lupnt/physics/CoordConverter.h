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

#include "Cheby.h"
#include "lupnt/core/Constants.h"

namespace ad = autodiff;

namespace LPT {

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

  static ad::VectorXreal Convert(const ad::VectorXreal rv_in,
                                 const ad::real epoch,
                                 const std::string coord_sys_in,
                                 const std::string coord_sys_out);
  static ad::VectorXreal Convert(const ad::VectorXreal rv_in,
                                 const ad::real epoch,
                                 const CoordSystem coord_sys_in,
                                 const CoordSystem coord_sys_out);

 private:
  static ad::MatrixXreal ComputeITRFtoGCRF(const ad::real tai);
  static ad::Matrix3real R1(ad::real phi);
  static ad::Matrix3real R2(ad::real phi);
  static ad::Matrix3real R3(ad::real phi);
  static ad::Matrix3real Skew(ad::Vector3real vec);
};
}  // namespace LPT