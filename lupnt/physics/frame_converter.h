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
  MI,           // Moon-centered Inertial Frame (Axis aligened with ICRF)
  PA,           // Moon-Fixed with principal axes
  ME,           // Moon-Fixed with mean-Earth / polar axes
  OP,           // Earth Orbit Frame
};

class FrameConverter {
 public:
  // Vector = func(real, Vector)
  static Vector6 Convert(real epoch, Vector6 rv_in, Frame frame_in,
                         Frame frame_out);
  static Vector3 Convert(real epoch, Vector3 r_in, Frame frame_in,
                         Frame frame_out);
  // Matrix = func(real, Matrix)
  static Matrix<-1, 6> Convert(real epoch, const Matrix<-1, 6> &rv_in,
                               Frame frame_in, Frame frame_out);
  static Matrix<-1, 3> Convert(real epoch, const Matrix<-1, 3> &r_in,
                               Frame frame_in, Frame frame_out);
  // Matrix = func(Vector, Vector)
  static Matrix<-1, 6> Convert(VectorX epoch, const Vector6 &rv_in,
                               Frame frame_in, Frame frame_out);
  static Matrix<-1, 3> Convert(VectorX epoch, const Vector3 &r_in,
                               Frame frame_in, Frame frame_out);
  // Matrix = func(Vector, Matrix)
  static Matrix<-1, 6> Convert(VectorX epoch, const Matrix<-1, 6> &rv_in,
                               Frame frame_in, Frame frame_out);
  static Matrix<-1, 3> Convert(VectorX epoch, const Matrix<-1, 3> &r_in,
                               Frame frame_in, Frame frame_out);

  static CartesianOrbitState Convert(real epoch,
                                     const CartesianOrbitState &state_in,
                                     Frame frame_out);

 private:
  static Matrix6 ComputeOpToMi(real epoch);
};
}  // namespace lupnt