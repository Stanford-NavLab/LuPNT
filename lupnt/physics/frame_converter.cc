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

#include "frame_converter.h"

#include <filesystem>

#include "cheby.h"
#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/orbit_state.h"
#include "spice_interface.h"

namespace lupnt {
CartesianOrbitState ConvertFrame(Real tai, const CartesianOrbitState &state_in,
                                 Frame frame_out) {
  Vec6 rv_in = state_in.GetVec();
  Vec6 rv_out = ConvertFrame(tai, rv_in, state_in.GetCoordSystem(), frame_out);
  return CartesianOrbitState(rv_out, frame_out);
}

Vec3 ConvertFrame(Real tai, const Vec3 &rv_in, Frame frame_in,
                  Frame frame_out) {
  Vec6 rv_in_6;
  rv_in_6 << rv_in, Vec3::Zero();
  Vec6 rv_out_6 = ConvertFrame(tai, rv_in_6, frame_in, frame_out);
  return rv_out_6.head(3);
}

Mat<-1, 6> ConvertFrame(Real tai, const Mat<-1, 6> &rv_in, Frame frame_in,
                        Frame frame_out) {
  Mat<-1, 6> rv_out(rv_in.rows(), 6);
  for (int i = 0; i < rv_in.rows(); i++) {
    rv_out.row(i) =
        ConvertFrame(tai, rv_in.row(i).transpose().eval(), frame_in, frame_out);
  }
  return rv_out;
}

Mat<-1, 3> ConvertFrame(Real tai, const Mat<-1, 3> &r_in, Frame frame_in,
                        Frame frame_out) {
  Mat<-1, 6> rv_in(r_in.rows(), 6);
  rv_in << r_in, Mat<-1, 3>::Zero(r_in.rows(), 3);
  Mat<-1, 6> rv_out = ConvertFrame(tai, rv_in, frame_in, frame_out);
  return rv_out.leftCols(3);
}

Mat<-1, 6> ConvertFrame(VecX tai, const Vec6 &rv_in, Frame frame_in,
                        Frame frame_out) {
  Mat<-1, 6> rv_out(tai.size(), 6);
  for (int i = 0; i < tai.size(); i++) {
    rv_out.row(i) = ConvertFrame(tai(i), rv_in, frame_in, frame_out);
  }
  return rv_out;
}

Mat<-1, 3> ConvertFrame(VecX tai, const Vec3 &r_in, Frame frame_in,
                        Frame frame_out) {
  Mat<-1, 6> rv_in(tai.size(), 6);
  rv_in << r_in, Mat<-1, 3>::Zero(tai.size(), 3);
  Mat<-1, 6> rv_out = ConvertFrame(tai, rv_in, frame_in, frame_out);
  return rv_out.leftCols(3);
}

Mat<-1, 6> ConvertFrame(VecX tai, const Mat<-1, 6> &rv_in, Frame frame_in,
                        Frame frame_out) {
  assert(tai.size() == rv_in.rows() && "Epoch and rv_in must have same size");
  Mat<-1, 6> rv_out(tai.size(), 6);
  for (int i = 0; i < tai.size(); i++) {
    rv_out.row(i) = ConvertFrame(tai(i), rv_in.row(i).transpose().eval(),
                                 frame_in, frame_out);
  }
  return rv_out;
}

Mat<-1, 3> ConvertFrame(VecX tai, const Mat<-1, 3> &r_in, Frame frame_in,
                        Frame frame_out) {
  assert(tai.size() == r_in.rows() && "Epoch and r_in must have same size");
  Mat<-1, 6> rv_in(tai.size(), 6);
  rv_in << r_in, Mat<-1, 3>::Zero(tai.size(), 3);
  Mat<-1, 6> rv_out = ConvertFrame(tai, rv_in, frame_in, frame_out);
  return rv_out.leftCols(3);
}

///
/// @brief Convert the state vector from one coordinate system to another (with
/// integer ID input)
///
/// @param tai Epoch of the state vector
/// @param rv_in   State vector in the original coordinate system
/// @param frame_in Coordinate system of the original state vector
/// @param frame_out Coordinate system of the converted state vector
/// @return Vec6  State vector in the converted coordinate system
///
Vec6 ConvertFrame(Real tai, const Vec6 &rv_in, Frame frame_in,
                  Frame frame_out) {
  if (frame_in == frame_out) {
    return rv_in;
  }

  switch (frame_in) {
    case ITRF: {
      switch (frame_out) {
        case GCRF: {
          // Convert to GCRF
          Mat6 Rrv_itrf2gcrf =
              GetFrameConversionMat(tai, Frame::ITRF, Frame::GCRF);
          Vec6 rv_gcrf = Rrv_itrf2gcrf * rv_in;
          return rv_gcrf;
        }
        default: {  // First convert to GCRF and then to the desired frame
          Vec6 rv_gcrf = ConvertFrame(tai, rv_in, ITRF, GCRF);
          Vec6 rv_out = ConvertFrame(tai, rv_gcrf, GCRF, frame_out);
          return rv_out;
        }
      }
    }

    case MOON_ME: {
      switch (frame_out) {
        case MOON_PA: {  // Convert to MOON_PA
          Vec3 r_ME = rv_in.head(3);
          Vec3 v_ME = rv_in.tail(3);

          // Rotation Mat MOON_ME (in DE421) -> MOON_PA (in DE440)
          // Reference:
          // https://iopscience.iop.org/article/10.3847/1538-3881/abd414/pdf
          Mat3 B_M = RotX(-0.2785 * DEG_PER_ARCSEC) *
                     RotY(-78.6944 * DEG_PER_ARCSEC) *
                     RotZ(-67.8526 * DEG_PER_ARCSEC);
          Mat3 B_M_inv = B_M.transpose();

          Vec3 r_PA = B_M * r_ME;
          Vec3 v_PA = B_M * v_ME;
          Vec6 rv_PA;
          rv_PA << r_PA, v_PA;
          return rv_PA;
        }
        default: {  // first convert to MOON_PA and then to the desired frame
          Vec6 rv_pa = ConvertFrame(tai, rv_in, MOON_ME, MOON_PA);
          Vec6 rv_out = ConvertFrame(tai, rv_pa, MOON_PA, frame_out);
          return rv_out;
        }
      }
    }

    case MOON_PA: {
      switch (frame_out) {
        case MOON_ME: {
          // Rotation Mat MOON_ME (in DE421) -> MOON_PA (in DE440)
          // Reference:
          // https://iopscience.iop.org/article/10.3847/1538-3881/abd414/pdf
          Vec3 r_PA = rv_in.head(3);
          Vec3 v_PA = rv_in.tail(3);
          Mat3 B_M = RotX(-0.2785 * DEG_PER_ARCSEC) *
                     RotY(-78.6944 * DEG_PER_ARCSEC) *
                     RotZ(-67.8526 * DEG_PER_ARCSEC);
          Vec3 r_ME = B_M * r_PA;
          Vec3 v_ME = B_M * v_PA;
          Vec6 rv_me;
          rv_me << r_ME, v_ME;
          return rv_me;
        }
        case MOON_CI: {  // Convert to Moon Inertial
          Mat6 Mrot = GetFrameConversionMat(tai, Frame::MOON_PA, Frame::GCRF);
          Vec6 rv_mi = Mrot * rv_in;
          return rv_mi;
        }
        default: {  // first convert to MOON_CI and then to the desired frame
          Vec6 rv_mi = ConvertFrame(tai, rv_in, MOON_PA, MOON_CI);
          Vec6 rv_out = ConvertFrame(tai, rv_mi, MOON_CI, frame_out);
          return rv_out;
        }
      }
    }

    case GCRF: {
      switch (frame_out) {
        case ICRF: {
          Vec6 rv_icrf_ssb2e = GetBodyPosVel(
              tai, NaifId::SOLAR_SYSTEM_BARYCENTER, NaifId::EARTH, Frame::GCRF);
          Vec6 rv_icrf = rv_in + rv_icrf_ssb2e;
          return rv_icrf;
        }
        case ITRF: {
          Mat6 Rrv_gcrf2itrf =
              GetFrameConversionMat(tai, Frame::GCRF, Frame::ITRF);
          Vec6 rv_itrf = Rrv_gcrf2itrf * rv_in;
          return rv_itrf;
        }
        case MOON_CI: {
          Vec6 rv_icrf_m2e =
              GetBodyPosVel(tai, NaifId::MOON, NaifId::EARTH, Frame::GCRF);
          Vec6 rv_mi = rv_in + rv_icrf_m2e;
          return rv_mi;
        }
        case EMR: {
          Vec6 rv_icrf_emb2e = GetBodyPosVel(tai, NaifId::EARTH_MOON_BARYCENTER,
                                             NaifId::EARTH, Frame::GCRF);
          Vec6 rv_emr = Inertial2Rtn(rv_icrf_emb2e, rv_in);
          return rv_emr;
        }
        case MOON_PA:
        case MOON_ME: {
          Vec6 rv_mi = ConvertFrame(tai, rv_in, GCRF, MOON_CI);
          Vec6 rv_out = ConvertFrame(tai, rv_mi, MOON_CI, frame_out);
          return rv_out;
        }
        default:
          assert(false && "Conversion not found");
      }
    }

    case MOON_CI: {
      switch (frame_out) {
        case GCRF: {
          Vec6 rv_icrf_e2m =
              GetBodyPosVel(tai, NaifId::EARTH, NaifId::MOON, Frame::MOON_CI);
          Vec6 rv_gcrf = rv_in + rv_icrf_e2m;
          return rv_gcrf;
        }
        case MOON_PA: {
          Mat6 Mrot = GetFrameConversionMat(tai, Frame::GCRF, Frame::MOON_PA);
          Vec6 rv_pa = Mrot * rv_in;
          return rv_pa;
        }
        case MOON_ME: {
          Vec6 rv_pa = ConvertFrame(tai, rv_in, MOON_CI, MOON_PA);
          Vec6 rv_out = ConvertFrame(tai, rv_pa, MOON_PA, frame_out);
          return rv_out;
        }
        case MOON_OP: {
          Mat6 R_op2mi = Op2Mi(tai);
          Vec6 rv_op = R_op2mi.transpose() * rv_in;
          return rv_op;
        }
        default: {
          Vec6 rv_gcrf = ConvertFrame(tai, rv_in, MOON_CI, GCRF);
          Vec6 rv_out = ConvertFrame(tai, rv_gcrf, GCRF, frame_out);
          return rv_out;
        }
      }
    }

    case ICRF: {
      switch (frame_out) {
        case GCRF: {
          Vec6 rv_icrf_ssb2e = GetBodyPosVel(
              tai, NaifId::SOLAR_SYSTEM_BARYCENTER, NaifId::EARTH, Frame::ICRF);
          Vec6 rv_gcrf = rv_in - rv_icrf_ssb2e;
          return rv_gcrf;
        }
        default: {
          Vec6 rv_gcrf = ConvertFrame(tai, rv_in, ICRF, GCRF);
          Vec6 rv_out = ConvertFrame(tai, rv_gcrf, GCRF, frame_out);
          return rv_out;
        }
      }
    }

    case EMR: {
      switch (frame_out) {
        case GCRF: {
          Vec6 rv_icrf_emb2e = GetBodyPosVel(
              tai, NaifId::EARTH, NaifId::EARTH_MOON_BARYCENTER, Frame::GCRF);
          Vec6 rv_gcrf = Rtn2Inertial(rv_icrf_emb2e, rv_in);
          return rv_gcrf;
        }
        default: {
          Vec6 rv_gcrf = ConvertFrame(tai, rv_in, EMR, GCRF);
          Vec6 rv_out = ConvertFrame(tai, rv_gcrf, GCRF, frame_out);
          return rv_out;
        }
      }
    }

    case MOON_OP: {
      switch (frame_out) {
        case MOON_CI: {
          Mat6 R_op2mi = Op2Mi(tai);
          Vec6 rv_mi = R_op2mi * rv_in;
          return rv_mi;
        }
        default: {
          Vec6 rv_mi = ConvertFrame(tai, rv_in, MOON_OP, MOON_CI);
          Vec6 rv_out = ConvertFrame(tai, rv_mi, MOON_CI, frame_out);
          return rv_out;
        }
      }
    }

    default: {
      assert(false && "Conversion not found");
    }
  }
}

Mat6 Op2Mi(Real tai) {
  Vec6 rv_earth_icrf =
      GetBodyPosVel(tai, NaifId::SUN, NaifId::EARTH, Frame::ICRF);
  Vec6 rv_moon_icrf =
      GetBodyPosVel(tai, NaifId::SUN, NaifId::MOON, Frame::ICRF);

  // Moon axis
  MatX iau_moon2icrf;
  iau_moon2icrf =
      GetFrameConversionMat(tai, Frame::MOON_PA, Frame::GCRF).block(0, 0, 3, 3);

  // IAU pole
  Vec3 iau_pole{0, 0, 1};
  Vec3 iau_pole_icrf = iau_moon2icrf * iau_pole;

  // MOON_OP unit vectors
  Vec3 dr = rv_earth_icrf.head(3) - rv_moon_icrf.head(3);
  Vec3 dv = rv_earth_icrf.tail(3) - rv_moon_icrf.tail(3);
  Vec3 z_op = dr.cross(dv).normalized();
  Vec3 x_op = iau_pole_icrf.cross(z_op).normalized();
  Vec3 y_op = z_op.cross(x_op).normalized();

  // create rotation matrix from MOON_OP to MOON_CI
  Mat3 R_op2mi;
  R_op2mi << x_op, y_op, z_op;

  Mat6 R_op2mi_tot = Mat6::Zero();
  R_op2mi_tot.topLeftCorner(3, 3) = R_op2mi;
  R_op2mi_tot.bottomRightCorner(3, 3) = R_op2mi;

  return R_op2mi_tot.cast<double>();
}

}  // namespace lupnt