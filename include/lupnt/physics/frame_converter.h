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

#include <functional>
#include <map>
#include <memory>

#include "cheby.h"
#include "lupnt/core/constants.h"

namespace lupnt {

  class CartesianOrbitState;

  enum Frame {
    ITRF,         // International Terrestrial Reference Frame
    ECEF = ITRF,  // Earth-Centered Earth-Fixed
    GCRF,         // Geocentric Reference System
    EME,          // Earth-Centered mean equator and equinox at J2000
    ECI = EME,    // Earth-Centered Inertial
    ICRF,         // International Celestial Reference System
    SER,          // Sun-Earth Rotating Frame
    GSE,          // Geocentric Solar Ecliptic
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

  extern std::map<std::pair<Frame, Frame>, std::function<Vec6(Real, const Vec6 &rv)>>
      frame_conversions;

  // Vec = func(real, Vec)
  Vec6 ConvertFrame(Real t_tai, const Vec6 &rv_in, Frame frame_in, Frame frame_out);
  Vec3 ConvertFrame(Real t_tai, const Vec3 &r_in, Frame frame_in, Frame frame_out);

  // Mat = func(real, Mat)
  Mat<-1, 6> ConvertFrame(Real t_tai, const Mat<-1, 6> &rv_in, Frame frame_in, Frame frame_out);
  Mat<-1, 3> ConvertFrame(Real t_tai, const Mat<-1, 3> &r_in, Frame frame_in, Frame frame_out);

  // Mat = func(Vec, Vec)
  Mat<-1, 6> ConvertFrame(VecX t_tai, const Vec6 &rv_in, Frame frame_in, Frame frame_out);
  Mat<-1, 3> ConvertFrame(VecX t_tai, const Vec3 &r_in, Frame frame_in, Frame frame_out);

  // Mat = func(Vec, Mat)
  Mat<-1, 6> ConvertFrame(VecX t_tai, const Mat<-1, 6> &rv_in, Frame frame_in, Frame frame_out);
  Mat<-1, 3> ConvertFrame(VecX t_tai, const Mat<-1, 3> &r_in, Frame frame_in, Frame frame_out);

  CartesianOrbitState ConvertFrame(Real t_tai, const CartesianOrbitState &state_in,
                                   Frame frame_out);

  Vec6 ITRF2GCRF(Real t_tai, const Vec6 &rv);
  Vec6 GCRF2ITRF(Real t_tai, const Vec6 &rv);

  Vec6 GCRF2EME(Real t_tai, const Vec6 &rv);
  Vec6 EME2GCRF(Real t_tai, const Vec6 &rv);

  Vec6 GCRF2ICRF(Real t_tai, const Vec6 &rv);
  Vec6 ICRF2GCRF(Real t_tai, const Vec6 &rv);

  Vec6 GCRF2MoonCI(Real t_tai, const Vec6 &rv);
  Vec6 MoonMI2GCRF(Real t_tai, const Vec6 &rv);

  Vec6 MoonMI2MoonPA(Real t_tai, const Vec6 &rv);
  Vec6 MoonPA2MoonCI(Real t_tai, const Vec6 &rv);

  Vec6 MoonPA2MoonME(Real t_tai, const Vec6 &rv);
  Vec6 MoonME2MoonPA(Real t_tai, const Vec6 &rv);

  Vec6 GCRF2EMR(Real t_tai, const Vec6 &rv);
  Vec6 EMR2GCRF(Real t_tai, const Vec6 &rv);

  Vec6 MoonME2MoonOP(Real t_tai, const Vec6 &rv);
  Vec6 MoonOP2MoonME(Real t_tai, const Vec6 &rv);

}  // namespace lupnt
