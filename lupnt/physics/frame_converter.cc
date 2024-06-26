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
#include "lupnt/physics/orbit_state_utils.h"
#include "spice_interface.h"

namespace lupnt {
CartesianOrbitState FrameConverter::Convert(real t_tai,
                                            const CartesianOrbitState &state_in,
                                            Frame frame_out) {
  Vector6 rv_in = state_in.GetVector();
  Vector6 rv_out = Convert(t_tai, rv_in, state_in.GetCoordSystem(), frame_out);
  return CartesianOrbitState(rv_out, frame_out);
}

Vector3 FrameConverter::Convert(real t_tai, Vector3 rv_in, Frame frame_in,
                                Frame frame_out) {
  Vector6 rv_in_6;
  rv_in_6 << rv_in, Vector3::Zero();
  Vector6 rv_out_6 = Convert(t_tai, rv_in_6, frame_in, frame_out);
  return rv_out_6.head(3);
}

Matrix<-1, 6> FrameConverter::Convert(real t_tai, const Matrix<-1, 6> &rv_in,
                                      Frame frame_in, Frame frame_out) {
  Matrix<-1, 6> rv_out(rv_in.rows(), 6);
  for (int i = 0; i < rv_in.rows(); i++) {
    rv_out.row(i) =
        Convert(t_tai, rv_in.row(i).transpose().eval(), frame_in, frame_out);
  }
  return rv_out;
}

Matrix<-1, 3> FrameConverter::Convert(real t_tai, const Matrix<-1, 3> &r_in,
                                      Frame frame_in, Frame frame_out) {
  Matrix<-1, 6> rv_in(r_in.rows(), 6);
  rv_in << r_in, Matrix<-1, 3>::Zero(r_in.rows(), 3);
  Matrix<-1, 6> rv_out = Convert(t_tai, rv_in, frame_in, frame_out);
  return rv_out.leftCols(3);
}

Matrix<-1, 6> FrameConverter::Convert(VectorX t_tai, const Vector6 &rv_in,
                                      Frame frame_in, Frame frame_out) {
  Matrix<-1, 6> rv_out(t_tai.size(), 6);
  for (int i = 0; i < t_tai.size(); i++) {
    rv_out.row(i) = Convert(t_tai(i), rv_in, frame_in, frame_out);
  }
  return rv_out;
}

Matrix<-1, 3> FrameConverter::Convert(VectorX t_tai, const Vector3 &r_in,
                                      Frame frame_in, Frame frame_out) {
  Matrix<-1, 6> rv_in(t_tai.size(), 6);
  rv_in << r_in, Matrix<-1, 3>::Zero(t_tai.size(), 3);
  Matrix<-1, 6> rv_out = Convert(t_tai, rv_in, frame_in, frame_out);
  return rv_out.leftCols(3);
}

Matrix<-1, 6> FrameConverter::Convert(VectorX t_tai, const Matrix<-1, 6> &rv_in,
                                      Frame frame_in, Frame frame_out) {
  assert(t_tai.size() == rv_in.rows() && "Epoch and rv_in must have same size");
  Matrix<-1, 6> rv_out(t_tai.size(), 6);
  for (int i = 0; i < t_tai.size(); i++) {
    rv_out.row(i) =
        Convert(t_tai(i), rv_in.row(i).transpose().eval(), frame_in, frame_out);
  }
  return rv_out;
}

Matrix<-1, 3> FrameConverter::Convert(VectorX t_tai, const Matrix<-1, 3> &r_in,
                                      Frame frame_in, Frame frame_out) {
  assert(t_tai.size() == r_in.rows() && "Epoch and r_in must have same size");
  Matrix<-1, 6> rv_in(t_tai.size(), 6);
  rv_in << r_in, Matrix<-1, 3>::Zero(t_tai.size(), 3);
  Matrix<-1, 6> rv_out = Convert(t_tai, rv_in, frame_in, frame_out);
  return rv_out.leftCols(3);
}

/**
 * @brief Convert the state vector from one coordinate system to another (with
 * integer ID input)
 *
 * @param t_tai Epoch of the state vector
 * @param rv_in   State vector in the original coordinate system
 * @param frame_in Coordinate system of the original state vector
 * @param frame_out Coordinate system of the converted state vector
 * @return Vector6  State vector in the converted coordinate system
 */
Vector6 FrameConverter::Convert(real t_tai, Vector6 rv_in, Frame frame_in,
                                Frame frame_out) {
  if (frame_in == frame_out) {
    return rv_in;
  }

  switch (frame_in) {
    case ITRF: {
      switch (frame_out) {
        case GCRF: {
          // Convert to GCRF
          Matrix6 Rrv_itrf2gcrf =
              GetFrameConversionMatrix(t_tai, Frame::ITRF, Frame::GCRF);
          Vector6 rv_gcrf = Rrv_itrf2gcrf * rv_in;
          return rv_gcrf;
        }
        default: {  // First convert to GCRF and then to the desired frame
          Vector6 rv_gcrf = Convert(t_tai, rv_in, ITRF, GCRF);
          Vector6 rv_out = Convert(t_tai, rv_gcrf, GCRF, frame_out);
          return rv_out;
        }
      }
    }

    case MOON_ME: {
      switch (frame_out) {
        case MOON_PA: {  // Convert to MOON_PA
          Vector3 r_ME = rv_in.head(3);
          Vector3 v_ME = rv_in.tail(3);

          // Rotation Matrix MOON_ME (in DE421) -> MOON_PA (in DE440)
          // Reference:
          // https://iopscience.iop.org/article/10.3847/1538-3881/abd414/pdf
          Matrix3 B_M = Rot1(-0.2785 * DEG_PER_ARCSEC) *
                        Rot2(-78.6944 * DEG_PER_ARCSEC) *
                        Rot3(-67.8526 * DEG_PER_ARCSEC);
          Matrix3 B_M_inv = B_M.transpose();

          Vector3 r_PA = B_M * r_ME;
          Vector3 v_PA = B_M * v_ME;
          Vector6 rv_PA;
          rv_PA << r_PA, v_PA;
          return rv_PA;
        }
        default: {  // first convert to MOON_PA and then to the desired frame
          Vector6 rv_pa = Convert(t_tai, rv_in, MOON_ME, MOON_PA);
          Vector6 rv_out = Convert(t_tai, rv_pa, MOON_PA, frame_out);
          return rv_out;
        }
      }
    }

    case MOON_PA: {
      switch (frame_out) {
        case MOON_ME: {
          // Rotation Matrix MOON_ME (in DE421) -> MOON_PA (in DE440)
          // Reference:
          // https://iopscience.iop.org/article/10.3847/1538-3881/abd414/pdf
          Vector3 r_PA = rv_in.head(3);
          Vector3 v_PA = rv_in.tail(3);
          Matrix3 B_M = Rot1(-0.2785 * DEG_PER_ARCSEC) *
                        Rot2(-78.6944 * DEG_PER_ARCSEC) *
                        Rot3(-67.8526 * DEG_PER_ARCSEC);
          Vector3 r_ME = B_M * r_PA;
          Vector3 v_ME = B_M * v_PA;
          Vector6 rv_me;
          rv_me << r_ME, v_ME;
          return rv_me;
        }
        case MOON_CI: {  // Convert to Moon Inertial
          Matrix6 Mrot =
              GetFrameConversionMatrix(t_tai, Frame::MOON_PA, Frame::GCRF);
          Vector6 rv_mi = Mrot * rv_in;
          return rv_mi;
        }
        default: {  // first convert to MOON_CI and then to the desired frame
          Vector6 rv_mi = Convert(t_tai, rv_in, MOON_PA, MOON_CI);
          Vector6 rv_out = Convert(t_tai, rv_mi, MOON_CI, frame_out);
          return rv_out;
        }
      }
    }

    case GCRF: {
      switch (frame_out) {
        case ICRF: {
          Vector6 rv_icrf_ssb2e =
              GetBodyPosVel(t_tai, NaifId::SOLAR_SYSTEM_BARYCENTER,
                            NaifId::EARTH, Frame::GCRF);
          Vector6 rv_icrf = rv_in + rv_icrf_ssb2e;
          return rv_icrf;
        }
        case ITRF: {
          Matrix6 Rrv_gcrf2itrf =
              GetFrameConversionMatrix(t_tai, Frame::GCRF, Frame::ITRF);
          Vector6 rv_itrf = Rrv_gcrf2itrf * rv_in;
          return rv_itrf;
        }
        case MOON_CI: {
          Vector6 rv_icrf_m2e =
              GetBodyPosVel(t_tai, NaifId::MOON, NaifId::EARTH, Frame::GCRF);
          Vector6 rv_mi = rv_in + rv_icrf_m2e;
          return rv_mi;
        }
        case EMR: {
          Vector6 rv_icrf_emb2e = GetBodyPosVel(
              t_tai, NaifId::EARTH_MOON_BARYCENTER, NaifId::EARTH, Frame::GCRF);
          Vector6 rv_emr = InertialToRtn(rv_icrf_emb2e, rv_in);
          return rv_emr;
        }
        case MOON_PA:
        case MOON_ME: {
          Vector6 rv_mi = Convert(t_tai, rv_in, GCRF, MOON_CI);
          Vector6 rv_out = Convert(t_tai, rv_mi, MOON_CI, frame_out);
          return rv_out;
        }
        default:
          assert(false && "Conversion not found");
      }
    }

    case MOON_CI: {
      switch (frame_out) {
        case GCRF: {
          Vector6 rv_icrf_e2m =
              GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::MOON, Frame::MOON_CI);
          Vector6 rv_gcrf = rv_in + rv_icrf_e2m;
          return rv_gcrf;
        }
        case MOON_PA: {
          Matrix6 Mrot =
              GetFrameConversionMatrix(t_tai, Frame::GCRF, Frame::MOON_PA);
          Vector6 rv_pa = Mrot * rv_in;
          return rv_pa;
        }
        case MOON_ME: {
          Vector6 rv_pa = Convert(t_tai, rv_in, MOON_CI, MOON_PA);
          Vector6 rv_out = Convert(t_tai, rv_pa, MOON_PA, frame_out);
          return rv_out;
        }
        case MOON_OP: {
          Matrix6 R_op2mi = ComputeOpToMi(t_tai);
          Vector6 rv_op = R_op2mi.transpose() * rv_in;
          return rv_op;
        }
        default: {
          Vector6 rv_gcrf = Convert(t_tai, rv_in, MOON_CI, GCRF);
          Vector6 rv_out = Convert(t_tai, rv_gcrf, GCRF, frame_out);
          return rv_out;
        }
      }
    }

    case ICRF: {
      switch (frame_out) {
        case GCRF: {
          Vector6 rv_icrf_ssb2e =
              GetBodyPosVel(t_tai, NaifId::SOLAR_SYSTEM_BARYCENTER,
                            NaifId::EARTH, Frame::ICRF);
          Vector6 rv_gcrf = rv_in - rv_icrf_ssb2e;
          return rv_gcrf;
        }
        default: {
          Vector6 rv_gcrf = Convert(t_tai, rv_in, ICRF, GCRF);
          Vector6 rv_out = Convert(t_tai, rv_gcrf, GCRF, frame_out);
          return rv_out;
        }
      }
    }

    case EMR: {
      switch (frame_out) {
        case GCRF: {
          Vector6 rv_icrf_emb2e = GetBodyPosVel(
              t_tai, NaifId::EARTH, NaifId::EARTH_MOON_BARYCENTER, Frame::GCRF);
          Vector6 rv_gcrf = RtnToInertial(rv_icrf_emb2e, rv_in);
          return rv_gcrf;
        }
        default: {
          Vector6 rv_gcrf = Convert(t_tai, rv_in, EMR, GCRF);
          Vector6 rv_out = Convert(t_tai, rv_gcrf, GCRF, frame_out);
          return rv_out;
        }
      }
    }

    case MOON_OP: {
      switch (frame_out) {
        case MOON_CI: {
          Matrix6 R_op2mi = ComputeOpToMi(t_tai);
          Vector6 rv_mi = R_op2mi * rv_in;
          return rv_mi;
        }
        default: {
          Vector6 rv_mi = Convert(t_tai, rv_in, MOON_OP, MOON_CI);
          Vector6 rv_out = Convert(t_tai, rv_mi, MOON_CI, frame_out);
          return rv_out;
        }
      }
    }

    default: {
      assert(false && "Conversion not found");
    }
  }
}

Matrix6 FrameConverter::ComputeOpToMi(real t_tai) {
  Vector6 rv_earth_icrf =
      GetBodyPosVel(t_tai, NaifId::SUN, NaifId::EARTH, Frame::ICRF);
  Vector6 rv_moon_icrf =
      GetBodyPosVel(t_tai, NaifId::SUN, NaifId::MOON, Frame::ICRF);

  // Moon axis
  MatrixX iau_moon2icrf;
  iau_moon2icrf = GetFrameConversionMatrix(t_tai, Frame::MOON_PA, Frame::GCRF)
                      .block(0, 0, 3, 3);

  // IAU pole
  Vector3 iau_pole{0, 0, 1};
  Vector3 iau_pole_icrf = iau_moon2icrf * iau_pole;

  // MOON_OP unit vectors
  Vector3 dr = rv_earth_icrf.head(3) - rv_moon_icrf.head(3);
  Vector3 dv = rv_earth_icrf.tail(3) - rv_moon_icrf.tail(3);
  Vector3 z_op = dr.cross(dv).normalized();
  Vector3 x_op = iau_pole_icrf.cross(z_op).normalized();
  Vector3 y_op = z_op.cross(x_op).normalized();

  // create rotation matrix from MOON_OP to MOON_CI
  Matrix3 R_op2mi;
  R_op2mi << x_op, y_op, z_op;

  Matrix6 R_op2mi_tot = Matrix6::Zero();
  R_op2mi_tot.topLeftCorner(3, 3) = R_op2mi;
  R_op2mi_tot.bottomRightCorner(3, 3) = R_op2mi;

  return R_op2mi_tot.cast<double>();
}

}  // namespace lupnt