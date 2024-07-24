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

extern std::map<std::pair<Frame, Frame>,
                std::function<Vec6(Real, const Vec6 &rv)>>
    frame_conversions;

// Vec = func(real, Vec)
Vec6 ConvertFrame(Real t_tai, const Vec6 &rv_in, Frame frame_in,
                  Frame frame_out);
Vec3 ConvertFrame(Real t_tai, const Vec3 &r_in, Frame frame_in,
                  Frame frame_out);

// Mat = func(real, Mat)
Mat<-1, 6> ConvertFrame(Real t_tai, const Mat<-1, 6> &rv_in, Frame frame_in,
                        Frame frame_out);
Mat<-1, 3> ConvertFrame(Real t_tai, const Mat<-1, 3> &r_in, Frame frame_in,
                        Frame frame_out);

// Mat = func(Vec, Vec)
Mat<-1, 6> ConvertFrame(VecX t_tai, const Vec6 &rv_in, Frame frame_in,
                        Frame frame_out);
Mat<-1, 3> ConvertFrame(VecX t_tai, const Vec3 &r_in, Frame frame_in,
                        Frame frame_out);

// Mat = func(Vec, Mat)
Mat<-1, 6> ConvertFrame(VecX t_tai, const Mat<-1, 6> &rv_in, Frame frame_in,
                        Frame frame_out);
Mat<-1, 3> ConvertFrame(VecX t_tai, const Mat<-1, 3> &r_in, Frame frame_in,
                        Frame frame_out);

CartesianOrbitState ConvertFrame(Real t_tai,
                                 const CartesianOrbitState &state_in,
                                 Frame frame_out);

Mat6 RotOp2Mi(Real t_tai);

Vec6 ITRFtoGCRF(Real t_tai, const Vec6 &rv);
Vec6 GCRFtoITRF(Real t_tai, const Vec6 &rv);

Vec6 GCRFtoEME(Real t_tai, const Vec6 &rv);
Vec6 EMEtoGCRF(Real t_tai, const Vec6 &rv);

Vec6 GCRFtoICRF(Real t_tai, const Vec6 &rv);
Vec6 ICRFtoGCRF(Real t_tai, const Vec6 &rv);

Vec6 GCRFtoMI(Real t_tai, const Vec6 &rv);
Vec6 MItoGCRF(Real t_tai, const Vec6 &rv);

Vec6 MItoPA(Real t_tai, const Vec6 &rv);
Vec6 PAtoMI(Real t_tai, const Vec6 &rv);

Vec6 PAtoME(Real t_tai, const Vec6 &rv);
Vec6 MEtoPA(Real t_tai, const Vec6 &rv);

Vec6 GCRFtoEMR(Real t_tai, const Vec6 &rv);
Vec6 EMRtoGCRF(Real t_tai, const Vec6 &rv);

}  // namespace lupnt