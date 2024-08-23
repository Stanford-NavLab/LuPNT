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

#include "lupnt/physics/frame_converter.h"

#include <filesystem>

#include "lupnt/core/constants.h"
#include "lupnt/data/eop.h"
#include "lupnt/data/iau_sofa.h"
#include "lupnt/data/kernels.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/cheby.h"
#include "lupnt/physics/frame_converter_spice.h"
#include "lupnt/physics/orbit_state.h"
#include "lupnt/physics/time_converter.h"

#define FRAME_CONVERSION(from, to, func) \
  {{Frame::from, Frame::to}, [](Real t, const Vec6& rv) -> Vec6 { return func(t, rv); }}

namespace lupnt {

  std::map<Frame, NaifId> frame_centers = {
      // Earth
      {Frame::ITRF, NaifId::EARTH},
      {Frame::ECEF, NaifId::EARTH},
      {Frame::GCRF, NaifId::EARTH},
      {Frame::EME, NaifId::EARTH},
      {Frame::ECI, NaifId::EARTH},
      {Frame::ICRF, NaifId::SOLAR_SYSTEM_BARYCENTER},
      // Moon
      {Frame::MOON_ME, NaifId::MOON},
      {Frame::MOON_CI, NaifId::MOON},
      {Frame::MOON_PA, NaifId::MOON},
      // Solar System
      {Frame::MARS_FIXED, NaifId::MARS},
      {Frame::VENUS_FIXED, NaifId::VENUS},
  };

  std::map<std::pair<Frame, Frame>, std::function<Vec6(Real, const Vec6& rv)>> frame_conversions = {
      FRAME_CONVERSION(GCRF, ITRF, GCRF2ITRF),
      FRAME_CONVERSION(ITRF, GCRF, ITRF2GCRF),
      FRAME_CONVERSION(GCRF, EME, GCRF2EME),
      FRAME_CONVERSION(EME, GCRF, EME2GCRF),
      FRAME_CONVERSION(GCRF, ICRF, GCRF2ICRF),
      FRAME_CONVERSION(ICRF, GCRF, ICRF2GCRF),
      FRAME_CONVERSION(GCRF, MOON_CI, GCRF2MoonCI),
      FRAME_CONVERSION(MOON_CI, GCRF, MoonMI2GCRF),
      FRAME_CONVERSION(MOON_CI, MOON_PA, MoonMI2MoonPA),
      FRAME_CONVERSION(MOON_PA, MOON_CI, MoonPA2MoonCI),
      FRAME_CONVERSION(MOON_PA, MOON_ME, MoonPA2MoonME),
      FRAME_CONVERSION(MOON_ME, MOON_PA, MoonME2MoonPA),
      FRAME_CONVERSION(GCRF, EMR, GCRF2EMR),
      FRAME_CONVERSION(EMR, GCRF, EMR2GCRF),
      FRAME_CONVERSION(MOON_CI, MOON_OP, MoonCI2MoonOP),
      FRAME_CONVERSION(MOON_OP, MOON_CI, MoonOP2MoonCI),
  };

  std::ostream& operator<<(std::ostream& os, Frame frame) {
    switch (frame) {
      case Frame::ITRF:
        os << "ITRF";
        break;
      case Frame::GCRF:
        os << "GCRF";
        break;
      case Frame::EME:
        os << "EME";
        break;
      case Frame::ICRF:
        os << "ICRF";
        break;
      case Frame::SER:
        os << "SER";
        break;
      case Frame::GSE:
        os << "GSE";
        break;
      case Frame::MOD:
        os << "MOD";
        break;
      case Frame::TOD:
        os << "TOD";
        break;
      case Frame::EMR:
        os << "EMR";
        break;
      case Frame::MOON_CI:
        os << "MOON_CI";
        break;
      case Frame::MOON_PA:
        os << "MOON_PA";
        break;
      case Frame::MOON_ME:
        os << "MOON_ME";
        break;
      case Frame::MOON_OP:
        os << "MOON_OP";
        break;
      case Frame::MARS_FIXED:
        os << "MARS_FIXED";
        break;
      case Frame::VENUS_FIXED:
        os << "VENUS_FIXED";
        break;
    }
    return os;
  }

  /// @brief Convert the state vector from one coordinate system to another
  /// (with integer ID input)
  /// @param t_tai Epoch of the state vector
  /// @param rv_in   State vector in the original coordinate system
  /// @param frame_in Coordinate system of the original state vector
  /// @param frame_out Coordinate system of the converted state vector
  /// @return Vec6  State vector in the converted coordinate system
  ///
  Vec6 ConvertFrameBase(Real t_tai, const Vec6& rv_in, Frame frame_in, Frame frame_out) {
    if (frame_in == frame_out) return rv_in;
    std::vector<Frame> path = FindShortestPath(frame_in, frame_out, frame_conversions);
    Vec6 rv_out = rv_in;
    for (size_t i = 0; i < path.size() - 1; i++) {
      std::function<Vec6(Real, const Vec6& rv)> f = frame_conversions[{path[i], path[i + 1]}];
      rv_out = f(t_tai, rv_out);
    }
    return rv_out;
  }

  CartesianOrbitState ConvertFrame(Real t_tai, const CartesianOrbitState& state_in, Frame frame_out,
                                   bool rotate_only) {
    Vec6 rv_in = state_in.GetVec();
    Vec6 rv_out = ConvertFrame(t_tai, rv_in, state_in.GetFrame(), frame_out, rotate_only);
    return CartesianOrbitState(rv_out, frame_out);
  }

  Vec6 ConvertFrame(Real t_tai, const Vec6& rv_in, Frame frame_in, Frame frame_out,
                    bool rotate_only) {
    if (rotate_only) {
      Vec6 zeros = Vec6::Zero();
      Vec3 translation = ConvertFrameBase(t_tai, zeros, frame_in, frame_out).head(3);

      MatX6 rv_out_tmp(3, 6);
      for (int i = 0; i < 3; i++) {
        Vec6 rv_in_tmp = Vec6::Zero();
        rv_in_tmp(i) = 1;
        rv_out_tmp.row(i) = ConvertFrameBase(t_tai, rv_in_tmp, frame_in, frame_out).transpose();
      }

      // Rot contains the basis vectors of frame_in (row-stacked) expressed in frame_out
      Mat3 Rot = rv_out_tmp.block(0, 0, 3, 3).rowwise() - translation.transpose();
      Vec3 r_out = Rot.transpose() * rv_in.head(3);
      Vec3 v_out = Rot.transpose() * rv_in.tail(3);
      Vec6 rv_out;
      rv_out << r_out, v_out;
      return rv_out;
    }
    Vec6 rv_out_6 = ConvertFrameBase(t_tai, rv_in, frame_in, frame_out);
    return rv_out_6.head(3);
  }

  Vec3 ConvertFrame(Real t_tai, const Vec3& r_in, Frame frame_in, Frame frame_out,
                    bool rotate_only) {
    Vec6 rv_in;
    rv_in << r_in, Vec3::Zero();
    Vec6 rv_out_6 = ConvertFrame(t_tai, rv_in, frame_in, frame_out, rotate_only);
    return rv_out_6.head(3);
  }

  MatX6 ConvertFrame(Real t_tai, const MatX6& rv_in, Frame frame_in, Frame frame_out,
                     bool rotate_only) {
    MatX6 rv_out(rv_in.rows(), 6);
    for (int i = 0; i < rv_in.rows(); i++) {
      rv_out.row(i)
          = ConvertFrame(t_tai, rv_in.row(i).transpose().eval(), frame_in, frame_out, rotate_only);
    }
    return rv_out;
  }

  MatX3 ConvertFrame(Real t_tai, const MatX3& r_in, Frame frame_in, Frame frame_out,
                     bool rotate_only) {
    MatX3 r_out(r_in.rows(), 3);
    for (int i = 0; i < r_in.rows(); i++) {
      Vec3 r_in_ = r_in.row(i).transpose();
      r_out.row(i) = ConvertFrame(t_tai, r_in_, frame_in, frame_out, rotate_only).head(3);
    }
    return r_out;
  }

  MatX6 ConvertFrame(VecX t_tai, const Vec6& rv_in, Frame frame_in, Frame frame_out,
                     bool rotate_only) {
    MatX6 rv_out(t_tai.size(), 6);
    for (int i = 0; i < t_tai.size(); i++) {
      rv_out.row(i) = ConvertFrame(t_tai(i), rv_in, frame_in, frame_out, rotate_only);
    }
    return rv_out;
  }

  MatX3 ConvertFrame(VecX t_tai, const Vec3& r_in, Frame frame_in, Frame frame_out,
                     bool rotate_only) {
    MatX3 r_out(t_tai.size(), 3);
    for (int i = 0; i < t_tai.size(); i++)
      r_out.row(i) = ConvertFrame(t_tai(i), r_in, frame_in, frame_out, rotate_only).head(3);
    return r_out;
  }

  MatX6 ConvertFrame(VecX t_tai, const MatX6& rv_in, Frame frame_in, Frame frame_out,
                     bool rotate_only) {
    assert(t_tai.size() == rv_in.rows() && "Epoch and rv_in must have same size");
    MatX6 rv_out(t_tai.size(), 6);
    for (int i = 0; i < t_tai.size(); i++) {
      rv_out.row(i) = ConvertFrame(t_tai(i), rv_in.row(i).transpose().eval(), frame_in, frame_out,
                                   rotate_only);
    }
    return rv_out;
  }

  MatX3 ConvertFrame(VecX t_tai, const MatX3& r_in, Frame frame_in, Frame frame_out,
                     bool rotate_only) {
    assert(t_tai.size() == r_in.rows() && "Epoch and r_in must have same size");
    MatX3 r_out(t_tai.size(), 3);
    for (int i = 0; i < t_tai.size(); i++) {
      Vec3 r_in_ = r_in.row(i).transpose();
      r_out.row(i) = ConvertFrame(t_tai(i), r_in_, frame_in, frame_out, rotate_only);
    }
    return r_out;
  }
}  // namespace lupnt
