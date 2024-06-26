/**
 * @file FrameConverter.cpp
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

class CartesianOrbitState;

enum Frame {
  ITRF,         // International Terrestrial Reference Frame
  ECEF = ITRF,  // Earth-Centered Earth-Fixed
  GCRF,         // Geocentric Reference System
  ECI = GCRF,   // Earth-Centered Inertial
  ICRF,         // International Celestial Reference System
  SER,          // Sun-Earth Rotating Frame
  GSE,          // Geocentric Solar Ecliptic
  EME,          // Earth-Centered mean equator and equinox at J2000 epoch
  MOD,          // Mean of date equatorial system
  TOD,          // True of date equatorial system
  EMR,          // Earth-Moon Rotating Frame
  MOON_CI,      // Moon-centered Inertial Frame (Axis aligened with ICRF)
  MOON_PA,      // Moon-Fixed with principal axes
  MOON_ME,      // Moon-Fixed with mean-Earth / polar axes
  MOON_OP,      // Earth Orbit Frame
  MARS_FIXED,   // Mars fixed frame
  VENUS_FIXED,  // Venus fixed frame
};

class FrameConverter {
 public:
  // Vec = func(real, Vec)
  static Vec6 Convert(real epoch, Vec6 rv_in, Frame frame_in, Frame frame_out);
  static Vec3 Convert(real epoch, Vec3 r_in, Frame frame_in, Frame frame_out);
  // Mat = func(real, Mat)
  static Mat<-1, 6> Convert(real epoch, const Mat<-1, 6> &rv_in, Frame frame_in,
                            Frame frame_out);
  static Mat<-1, 3> Convert(real epoch, const Mat<-1, 3> &r_in, Frame frame_in,
                            Frame frame_out);
  // Mat = func(Vec, Vec)
  static Mat<-1, 6> Convert(VecX epoch, const Vec6 &rv_in, Frame frame_in,
                            Frame frame_out);
  static Mat<-1, 3> Convert(VecX epoch, const Vec3 &r_in, Frame frame_in,
                            Frame frame_out);
  // Mat = func(Vec, Mat)
  static Mat<-1, 6> Convert(VecX epoch, const Mat<-1, 6> &rv_in, Frame frame_in,
                            Frame frame_out);
  static Mat<-1, 3> Convert(VecX epoch, const Mat<-1, 3> &r_in, Frame frame_in,
                            Frame frame_out);

  static CartesianOrbitState Convert(real epoch,
                                     const CartesianOrbitState &state_in,
                                     Frame frame_out);

 private:
  static Mat6 ComputeOpToMi(real epoch);
};
}  // namespace lupnt