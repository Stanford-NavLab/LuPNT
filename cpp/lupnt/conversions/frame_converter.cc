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

#include "lupnt/conversions/frame_converter.h"

#include "lupnt/conversions/frame_conversions.h"
#include "lupnt/conversions/frame_converter.h"
#include "lupnt/conversions/frame_converter_spice.h"
#include "lupnt/conversions/state_converter.h"
#include "lupnt/core/constants.h"
#include "lupnt/states/state.h"

#define FRAME_CONVERSION(from, to, func) \
  {{Frame::from, Frame::to}, [](Real t, const Vec6& rv) -> Vec6 { return func(t, rv); }}

namespace lupnt {

  const std::map<Frame, BodyId> frame_centers = {
      // Earth
      {Frame::ITRF, BodyId::EARTH},
      {Frame::ECEF, BodyId::EARTH},
      {Frame::GCRF, BodyId::EARTH},
      {Frame::EME, BodyId::EARTH},
      {Frame::ECI, BodyId::EARTH},
      {Frame::ICRF, BodyId::SOLAR_SYSTEM_BARYCENTER},
      // Moon
      {Frame::MOON_ME, BodyId::MOON},
      {Frame::MOON_CI, BodyId::MOON},
      {Frame::MOON_PA, BodyId::MOON},
      {Frame::MOON_OP, BodyId::MOON},
      // Solar System -- body-fixed frames
      {Frame::MERCURY_FIXED, BodyId::MERCURY},
      {Frame::VENUS_FIXED, BodyId::VENUS},
      {Frame::MARS_FIXED, BodyId::MARS},
      {Frame::JUPITER_FIXED, BodyId::JUPITER},
      {Frame::SATURN_FIXED, BodyId::SATURN},
      {Frame::URANUS_FIXED, BodyId::URANUS},
      {Frame::NEPTUNE_FIXED, BodyId::NEPTUNE},
      // Solar System -- planet-centered inertial frames
      {Frame::MERCURY_CI, BodyId::MERCURY},
      {Frame::VENUS_CI, BodyId::VENUS},
      {Frame::MARS_CI, BodyId::MARS},
      {Frame::JUPITER_CI, BodyId::JUPITER},
      {Frame::SATURN_CI, BodyId::SATURN},
      {Frame::URANUS_CI, BodyId::URANUS},
      {Frame::NEPTUNE_CI, BodyId::NEPTUNE},
  };

  BodyId GetFrameCenter(Frame frame) {
    auto it = frame_centers.find(frame);
    LUPNT_CHECK(it != frame_centers.end(),
                fmt::format("Frame {} not found in frame_centers", enum_name(frame)),
                "FrameConverter");
    return it->second;
  }

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
      case Frame::MERCURY_FIXED: os << "MERCURY_FIXED"; break;
      case Frame::VENUS_FIXED: os << "VENUS_FIXED"; break;
      case Frame::MARS_FIXED: os << "MARS_FIXED"; break;
      case Frame::JUPITER_FIXED: os << "JUPITER_FIXED"; break;
      case Frame::SATURN_FIXED: os << "SATURN_FIXED"; break;
      case Frame::URANUS_FIXED: os << "URANUS_FIXED"; break;
      case Frame::NEPTUNE_FIXED: os << "NEPTUNE_FIXED"; break;
      case Frame::MERCURY_CI: os << "MERCURY_CI"; break;
      case Frame::VENUS_CI: os << "VENUS_CI"; break;
      case Frame::MARS_CI: os << "MARS_CI"; break;
      case Frame::JUPITER_CI: os << "JUPITER_CI"; break;
      case Frame::SATURN_CI: os << "SATURN_CI"; break;
      case Frame::URANUS_CI: os << "URANUS_CI"; break;
      case Frame::NEPTUNE_CI: os << "NEPTUNE_CI"; break;
      default: throw std::runtime_error("Frame not implemented");
    }
    return os;
  }

  /// @brief Convert the state vector from one coordinate system to another
  /// (with integer ID input)
  /// @param t_tdb Epoch of the state vector
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
  ///
  /// The solar-system planets hang off the ICRF (SSB) hub: each <PLANET>_CI is
  /// ICRF-aligned and centered on the planet, and <PLANET>_FIXED co-rotates with
  /// it via the IAU orientation model (see frame_conversions.cc).

  /// @brief Convert to/from a planet <PLANET>_CI / <PLANET>_FIXED frame by routing
  /// through the SSB-centered ICRF hub. Handles the case where the *other* endpoint
  /// is itself a (possibly different) planet frame, an Earth/Moon frame, or ICRF.
  template <int N>
  Vec<N> ConvertPlanetFrame(Real t_tdb, const Vec<N>& rv_in, Frame frame_in, Frame frame_out) {
    // Fast path: both frames belong to the same planet -> a pure rotation about
    // the shared center. Avoids round-tripping through the SSB (~1e11 m), which
    // would otherwise lose precision on the small body-relative state.
    if (IsPlanetFrame(frame_in) && IsPlanetFrame(frame_out)
        && GetFrameCenter(frame_in) == GetFrameCenter(frame_out)) {
      BodyId body = GetFrameCenter(frame_in);
      bool in_fixed = IsPlanetFixedFrame(frame_in);
      bool out_fixed = IsPlanetFixedFrame(frame_out);
      if (in_fixed == out_fixed) return rv_in;  // CI<->CI or FIXED<->FIXED (same body)
      if (out_fixed) return BodyCiToFixed(t_tdb, rv_in, body);
      return BodyFixedToCi(t_tdb, rv_in, body);
    }
    // 1) Bring the input to the ICRF hub.
    Vec<N> rv_icrf;
    if (IsPlanetCiFrame(frame_in)) {
      rv_icrf = PlanetCiToIcrf(t_tdb, rv_in, GetFrameCenter(frame_in));
    } else if (IsPlanetFixedFrame(frame_in)) {
      Vec<N> rv_ci = BodyFixedToCi(t_tdb, rv_in, GetFrameCenter(frame_in));
      rv_icrf = PlanetCiToIcrf(t_tdb, rv_ci, GetFrameCenter(frame_in));
    } else {
      rv_icrf = ConvertFrame(t_tdb, rv_in, frame_in, Frame::ICRF);  // non-planet endpoint
    }
    // 2) Emit from the ICRF hub to the requested output.
    if (IsPlanetCiFrame(frame_out)) {
      return IcrfToPlanetCi(t_tdb, rv_icrf, GetFrameCenter(frame_out));
    } else if (IsPlanetFixedFrame(frame_out)) {
      Vec<N> rv_ci = IcrfToPlanetCi(t_tdb, rv_icrf, GetFrameCenter(frame_out));
      return BodyCiToFixed(t_tdb, rv_ci, GetFrameCenter(frame_out));
    }
    return ConvertFrame(t_tdb, rv_icrf, Frame::ICRF, frame_out);  // non-planet endpoint
  }

  template <int N>
  Vec<N> ConvertFrameBase(Real t_tdb, const Vec<N>& rv_in, Frame frame_in, Frame frame_out) {
    if (frame_in == frame_out) return rv_in;
    // Solar-system planet CI/FIXED frames (generic IAU orientation, ICRF hub).
    if (IsPlanetFrame(frame_in) || IsPlanetFrame(frame_out))
      return ConvertPlanetFrame(t_tdb, rv_in, frame_in, frame_out);
    switch (frame_in) {
      // Earth ***************
      case Frame::ICRF:
        return ConvertFrame(t_tdb, IcrfToGcrf(t_tdb, rv_in), Frame::GCRF, frame_out);
      case Frame::ITRF:  // Frame::ECEF:
        return ConvertFrame(t_tdb, ItrfToGcrf(t_tdb, rv_in), Frame::GCRF, frame_out);
      case Frame::GCRF:
        switch (frame_out) {
          case Frame::ICRF: return GcrfToIcrf(t_tdb, rv_in);
          case Frame::ITRF: return GcrfToItrf(t_tdb, rv_in);
          case Frame::SER:
          case Frame::GSE: throw std::runtime_error("Not implemented");
          case Frame::EMR: return GcrfToEmr(t_tdb, rv_in);
          case Frame::EME: return GcrfToEme(rv_in);
          case Frame::MOD:
          case Frame::TOD: throw std::runtime_error("Not implemented");
          case Frame::MOON_PA:
          case Frame::MOON_ME:
          case Frame::MOON_OP:
          case Frame::MOON_CI:
            return ConvertFrame(t_tdb, GcrfToMoonCi(t_tdb, rv_in), Frame::MOON_CI, frame_out);
          default: throw std::runtime_error("Not implemented");
        }
      case Frame::EME:  // Frame::ECI:
        return ConvertFrame(t_tdb, EmeToGcrf(rv_in), Frame::GCRF, frame_out);
      case Frame::EMR: return ConvertFrame(t_tdb, EmrToGcrf(t_tdb, rv_in), Frame::GCRF, frame_out);
      case Frame::SER:
      case Frame::GSE:
      case Frame::MOD:
      case Frame::TOD: throw std::runtime_error("Not implemented");
      // Moon ***************
      case Frame::MOON_CI:
        switch (frame_out) {
          case Frame::MOON_OP: return MoonCiToOp(t_tdb, rv_in);
          case Frame::MOON_PA: return MoonCiToPa(t_tdb, rv_in);
          case Frame::MOON_ME: return MoonPaToMe(MoonCiToPa(t_tdb, rv_in));
          default: return ConvertFrame(t_tdb, MoonCiToGcrf(t_tdb, rv_in), Frame::GCRF, frame_out);
        }
      case Frame::MOON_ME: return ConvertFrame(t_tdb, MoonMeToPa(rv_in), Frame::MOON_PA, frame_out);
      case Frame::MOON_PA: {
        if (frame_out == Frame::MOON_ME) return MoonPaToMe(rv_in);
        return ConvertFrame(t_tdb, MoonPaToCi(t_tdb, rv_in), Frame::MOON_CI, frame_out);
      }
      case Frame::MOON_OP:
        return ConvertFrame(t_tdb, MoonOpToCi(t_tdb, rv_in), Frame::MOON_CI, frame_out);
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
      case Frame::UNDEFINED: throw std::runtime_error("Not implemented");
    }
    throw std::runtime_error("Conversion not implemented");
    return Vec<N>::Zero();
  }

  template Vec6 ConvertFrameBase(Real t_tdb, const Vec6& rv_in, Frame frame_in, Frame frame_out);
  template Vec3 ConvertFrameBase(Real t_tdb, const Vec3& r_in, Frame frame_in, Frame frame_out);

  Vec3 ConvertFrame(Real t_tdb, const Vec3& r_in, Frame frame_in, Frame frame_out,
                    bool rotate_only) {
    if (frame_in == frame_out) return r_in;
    if (!rotate_only) return ConvertFrameBase(t_tdb, r_in, frame_in, frame_out);
    auto [R, r] = GetFrameRotationTranslation(t_tdb, frame_in, frame_out);
    return R * r_in + r;
  }

  Vec6 ConvertFrame(Real t_tdb, const Vec6& rv_in, Frame frame_in, Frame frame_out,
                    bool rotate_only) {
    if (frame_in == frame_out) return rv_in;
    if (!rotate_only) return ConvertFrameBase(t_tdb, rv_in, frame_in, frame_out);
    auto [R, r] = GetFrameRotationTranslation(t_tdb, frame_in, frame_out);
    Vec6 rv_out;
    rv_out.head(3) = R * rv_in.head(3) + r;
    rv_out.tail(3) = R * rv_in.tail(3);
    return rv_out;
  }

  MatX6 ConvertFrame(Real t_tdb, const MatX6& rv_in, Frame frame_in, Frame frame_out,
                     bool rotate_only) {
    MatX6 rv_out(rv_in.rows(), 6);
    if (rotate_only) {
      auto [R, r] = GetFrameRotationTranslation(t_tdb, frame_in, frame_out);
      for (int i = 0; i < rv_in.rows(); i++) {
        rv_out.row(i).head(3) = R * rv_in.row(i).head(3).transpose() + r;
        rv_out.row(i).tail(3) = R * rv_in.row(i).tail(3).transpose();
      }
    } else {
      for (int i = 0; i < rv_in.rows(); i++)
        rv_out.row(i)
            = ConvertFrameBase(t_tdb, rv_in.row(i).transpose().eval(), frame_in, frame_out);
    }
    return rv_out;
  }

  MatX3 ConvertFrame(Real t_tdb, const MatX3& r_in, Frame frame_in, Frame frame_out,
                     bool rotate_only) {
    MatX3 r_out(r_in.rows(), 3);
    if (rotate_only) {
      auto [R, r] = GetFrameRotationTranslation(t_tdb, frame_in, frame_out);
      for (int i = 0; i < r_in.rows(); i++) r_out.row(i) = R * r_in.row(i).transpose() + r;
    } else {
      for (int i = 0; i < r_in.rows(); i++)
        r_out.row(i) = ConvertFrameBase(t_tdb, r_in.row(i).transpose().eval(), frame_in, frame_out);
    }
    return r_out;
  }

  MatX6 ConvertFrame(const VecX& t_tdb, const Vec6& rv_in, Frame frame_in, Frame frame_out,
                     bool rotate_only) {
    MatX6 rv_out(t_tdb.size(), 6);
    for (int i = 0; i < t_tdb.size(); i++)
      rv_out.row(i) = ConvertFrame(t_tdb(i), rv_in, frame_in, frame_out, rotate_only);
    return rv_out;
  }

  MatX3 ConvertFrame(const VecX& t_tdb, const Vec3& r_in, Frame frame_in, Frame frame_out,
                     bool rotate_only) {
    MatX3 r_out(t_tdb.size(), 3);
    for (int i = 0; i < t_tdb.size(); i++)
      r_out.row(i) = ConvertFrame(t_tdb(i), r_in, frame_in, frame_out, rotate_only);
    return r_out;
  }

  MatX6 ConvertFrame(const VecX& t_tdb, const MatX6& rv_in, Frame frame_in, Frame frame_out,
                     bool rotate_only) {
    LUPNT_CHECK(t_tdb.size() == rv_in.rows(), "Epoch and rv_in must have same size",
                "FrameConverter");
    MatX6 rv_out(t_tdb.size(), 6);
    for (int i = 0; i < t_tdb.size(); i++)
      rv_out.row(i) = ConvertFrame(t_tdb(i), rv_in.row(i).transpose().eval(), frame_in, frame_out,
                                   rotate_only);
    return rv_out;
  }

  MatX3 ConvertFrame(const VecX& t_tdb, const MatX3& r_in, Frame frame_in, Frame frame_out,
                     bool rotate_only) {
    if (t_tdb.size() != r_in.rows()) throw std::runtime_error("Epoch and r_in must have same size");
    MatX3 r_out(t_tdb.size(), 3);
    for (int i = 0; i < t_tdb.size(); i++) {
      Vec3 r_in_ = r_in.row(i).transpose();
      r_out.row(i) = ConvertFrame(t_tdb(i), r_in_, frame_in, frame_out, rotate_only);
    }
    return r_out;
  }

  State ConvertFrame(Real t_tdb, const State& x_in, Frame frame_out, bool rotate_only) {
    if (x_in.GetFrame() == frame_out) return x_in;
    Cart6 rv_in = ConvertState(x_in, Cart6::TYPE);
    Cart6 rv_out = ConvertFrame(t_tdb, Vec6(rv_in), x_in.GetFrame(), frame_out, rotate_only);
    rv_out.SetFrame(frame_out);
    return ConvertState(rv_out, x_in.GetType());
  }

  /// @brief r_to = R * r_from + r
  /// @param t_tdb
  /// @param from_frame
  /// @param to_frame
  /// @return
  std::pair<Mat3, Vec3> GetFrameRotationTranslation(Real t_tdb, Frame from_frame, Frame to_frame) {
    MatX3 I3 = Mat3::Identity();
    Mat3 M = ConvertFrame(t_tdb, I3, from_frame, to_frame);
    MatX3 r0 = MatX3::Zero(1, 3);
    Vec3 r = ConvertFrame(t_tdb, r0, from_frame, to_frame).transpose();
    Mat3 R = (M.rowwise() - r.transpose()).transpose();
    return std::make_pair(R, r);
  }
  /// @brief r_to = R * r_from + r
  /// @param t_tdb
  /// @param from_frame
  /// @param to_frame
  /// @return
  /// @note Returns a pair of rotation matrix and translation vector
  ///       The rotation matrix is the rotation from from_frame to to_frame
  ///       The translation vector is the translation from from_frame to to_frame
  /// Returns (R6, t6) such that x_to = R6 * x_from + t6
  /// R6 is 6x6; t6 is 6x1.
  std::pair<Mat6, Vec6> GetFrameRotationTranslationRv(Real t_tdb, Frame from_frame,
                                                      Frame to_frame) {
    // Transform the zero state -> this is the affine offset t6.
    MatX6 zero_row = MatX6::Zero(1, 6);
    Vec6 t6 = ConvertFrame(t_tdb, zero_row, from_frame, to_frame, false).row(0).transpose();  // 6x1

    // Transform the six basis states (each row is one basis vector).
    MatX6 I6 = Mat6::Identity();
    MatX6 Y = ConvertFrame(t_tdb, I6, from_frame, to_frame, false);  // 6x6

    // Build R6: columns are (Y_row_i^T - t6), for i=0..5.
    Mat6 R6;
    for (int j = 0; j < 6; ++j) {
      R6.col(j) = Y.row(j).transpose() - t6;  // 6x1
    }

    return {R6, t6};
  }

}  // namespace lupnt
