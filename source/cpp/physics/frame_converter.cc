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
#include "lupnt/physics/frame_conversions.h"
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

  std::ostream& operator<<(std::ostream& os, Frame frame) {
    switch (frame) {
      case Frame::ITRF: os << "ITRF"; break;
      case Frame::GCRF: os << "GCRF"; break;
      case Frame::EME: os << "EME"; break;
      case Frame::ICRF: os << "ICRF"; break;
      case Frame::SER: os << "SER"; break;
      case Frame::GSE: os << "GSE"; break;
      case Frame::MOD: os << "MOD"; break;
      case Frame::TOD: os << "TOD"; break;
      case Frame::EMR: os << "EMR"; break;
      case Frame::MOON_CI: os << "MOON_CI"; break;
      case Frame::MOON_PA: os << "MOON_PA"; break;
      case Frame::MOON_ME: os << "MOON_ME"; break;
      case Frame::MOON_OP: os << "MOON_OP"; break;
      case Frame::MARS_FIXED: os << "MARS_FIXED"; break;
      case Frame::VENUS_FIXED: os << "VENUS_FIXED"; break;
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
  /// @note
  ///     ITRF  TOD
  ///      |     |
  ///     TIRS  MOD   ME
  ///      |     |    |
  ///     CIRS  EME   PA
  ///        \  /     |
  /// ICRF -- GCRF -- MI -- OP
  ///       /  |  \  /
  ///     SER GSE EMR
  Vec6 ConvertFrameBase(Real t_tai, const Vec6& rv_in, Frame frame_in, Frame frame_out) {
    if (frame_in == frame_out) return rv_in;
    switch (frame_in) {
      // Earth ***************
      case Frame::ICRF:
        return ConvertFrameBase(t_tai, ICRF2GCRF(t_tai, rv_in), Frame::GCRF, frame_out);
      case Frame::ITRF:  // Frame::ECEF:
        return ConvertFrameBase(t_tai, ITRF2GCRF(t_tai, rv_in), Frame::GCRF, frame_out);
      case Frame::GCRF:
        switch (frame_out) {
          case Frame::ICRF: return GCRF2ICRF(t_tai, rv_in);
          case Frame::ITRF: return GCRF2ITRF(t_tai, rv_in);
          case Frame::SER:
          case Frame::GSE: throw std::runtime_error("Not implemented");
          case Frame::EMR: return GCRF2EMR(t_tai, rv_in);
          case Frame::EME: return GCRF2EME(t_tai, rv_in);
          case Frame::MOD:
          case Frame::TOD: throw std::runtime_error("Not implemented");
          case Frame::MOON_PA:
          case Frame::MOON_ME:
          case Frame::MOON_OP:
          case Frame::MOON_CI:
            return ConvertFrameBase(t_tai, GCRF2MoonCI(t_tai, rv_in), Frame::MOON_CI, frame_out);
          default: break;
        }
      case Frame::EME:  // Frame::ECI:
        return ConvertFrameBase(t_tai, EME2GCRF(t_tai, rv_in), Frame::GCRF, frame_out);
      case Frame::EMR:
        return ConvertFrameBase(t_tai, EMR2GCRF(t_tai, rv_in), Frame::GCRF, frame_out);
      case Frame::SER:
      case Frame::GSE:
      case Frame::MOD:
      case Frame::TOD: throw std::runtime_error("Not implemented");
      // Moon ***************
      case Frame::MOON_CI:
        switch (frame_out) {
          case Frame::MOON_OP: return MoonCI2MoonOP(t_tai, rv_in);
          case Frame::MOON_PA: return MoonCI2MoonPA(t_tai, rv_in);
          case Frame::MOON_ME: return MoonPA2MoonME(t_tai, MoonCI2MoonPA(t_tai, rv_in));
          default:
            return ConvertFrameBase(t_tai, MoonCI2GCRF(t_tai, rv_in), Frame::GCRF, frame_out);
        }
      case Frame::MOON_ME:
        return ConvertFrameBase(t_tai, MoonME2MoonPA(t_tai, rv_in), Frame::MOON_PA, frame_out);
      case Frame::MOON_PA: {
        if (frame_out == Frame::MOON_ME) return MoonPA2MoonME(t_tai, rv_in);
        return ConvertFrameBase(t_tai, MoonPA2MoonCI(t_tai, rv_in), Frame::MOON_CI, frame_out);
      }
      case Frame::MOON_OP:
        return ConvertFrameBase(t_tai, MoonOP2MoonCI(t_tai, rv_in), Frame::MOON_CI, frame_out);
        // Solar System ***************
      case Frame::MERCURY_FIXED:
      case Frame::VENUS_FIXED:
      case Frame::MARS_FIXED:
      case Frame::JUPITER_FIXED:
      case Frame::SATURN_FIXED:
      case Frame::URANUS_FIXED:
      case Frame::NEPTUNE_FIXED:
      case Frame::MERCURY_CI:
      case Frame::VENUS_CI:
      case Frame::MARS_CI:
      case Frame::JUPITER_CI:
      case Frame::SATURN_CI:
      case Frame::URANUS_CI:
      case Frame::NEPTUNE_CI:
      case Frame::NONE: throw std::runtime_error("Not implemented");
    }
    throw std::runtime_error("Conversion not implemented");
    return Vec6::Zero();
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

      // Rot contains the basis vectors of frame_in (row-stacked) expressed in
      // frame_out
      Mat3 Rot = rv_out_tmp.block(0, 0, 3, 3).rowwise() - translation.transpose();
      Vec3 r_out = Rot.transpose() * rv_in.head(3);
      Vec3 v_out = Rot.transpose() * rv_in.tail(3);
      Vec6 rv_out;
      rv_out << r_out, v_out;
      return rv_out;
    } else {
      Vec6 rv_out_6 = ConvertFrameBase(t_tai, rv_in, frame_in, frame_out);
      return rv_out_6;
    }
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

  Mat3 GetFrameConversionMatrix(Real t_tai, Frame from_frame, Frame to_frame) {
    MatX3 I = Mat3::Identity();
    return ConvertFrame(t_tai, I, to_frame, from_frame, true);
  }
}  // namespace lupnt
