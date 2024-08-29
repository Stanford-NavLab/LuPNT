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

#include "lupnt/physics/frame_converter_spice.h"

#include <filesystem>

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/cheby.h"
#include "lupnt/physics/orbit_state.h"
#include "lupnt/physics/spice_interface.h"

namespace lupnt {

  namespace spice {

    CartesianOrbitState ConvertFrameSpice(Real t_tai, const CartesianOrbitState& state_in,
                                          Frame frame_out) {
      Vec6 rv_in = state_in.GetVec();
      Vec6 rv_out = ConvertFrameSpice(t_tai, rv_in, state_in.GetFrame(), frame_out);
      return CartesianOrbitState(rv_out, frame_out);
    }

    Vec3 ConvertFrameSpice(Real t_tai, const Vec3& r_in, Frame frame_in, Frame frame_out) {
      Vec6 rv_in;
      rv_in << r_in, Vec3::Zero();
      Vec6 rv_out_6 = ConvertFrameSpice(t_tai, rv_in, frame_in, frame_out);
      return rv_out_6.head(3);
    }

    MatX6 ConvertFrameSpice(Real t_tai, const MatX6& rv_in, Frame frame_in, Frame frame_out) {
      MatX6 rv_out(rv_in.rows(), 6);
      for (int i = 0; i < rv_in.rows(); i++) {
        rv_out.row(i)
            = ConvertFrameSpice(t_tai, rv_in.row(i).transpose().eval(), frame_in, frame_out);
      }
      return rv_out;
    }

    MatX3 ConvertFrameSpice(Real t_tai, const MatX3& r_in, Frame frame_in, Frame frame_out) {
      MatX6 rv_in(r_in.rows(), 6);
      rv_in << r_in, MatX3::Zero(r_in.rows(), 3);
      MatX6 rv_out = ConvertFrameSpice(t_tai, rv_in, frame_in, frame_out);
      return rv_out.leftCols(3);
    }

    MatX6 ConvertFrameSpice(VecX t_tai, const Vec6& rv_in, Frame frame_in, Frame frame_out) {
      MatX6 rv_out(t_tai.size(), 6);
      for (int i = 0; i < t_tai.size(); i++) {
        rv_out.row(i) = ConvertFrameSpice(t_tai(i), rv_in, frame_in, frame_out);
      }
      return rv_out;
    }

    MatX3 ConvertFrameSpice(VecX t_tai, const Vec3& r_in, Frame frame_in, Frame frame_out) {
      MatX6 rv_in(t_tai.size(), 6);
      rv_in << r_in, MatX3::Zero(t_tai.size(), 3);
      MatX6 rv_out = ConvertFrameSpice(t_tai, rv_in, frame_in, frame_out);
      return rv_out.leftCols(3);
    }

    MatX6 ConvertFrameSpice(VecX t_tai, const MatX6& rv_in, Frame frame_in, Frame frame_out) {
      assert(t_tai.size() == rv_in.rows() && "Epoch and rv_in must have same size");
      MatX6 rv_out(t_tai.size(), 6);
      for (int i = 0; i < t_tai.size(); i++) {
        rv_out.row(i)
            = ConvertFrameSpice(t_tai(i), rv_in.row(i).transpose().eval(), frame_in, frame_out);
      }
      return rv_out;
    }

    MatX3 ConvertFrameSpice(VecX t_tai, const MatX3& r_in, Frame frame_in, Frame frame_out) {
      assert(t_tai.size() == r_in.rows() && "Epoch and r_in must have same size");
      MatX6 rv_in(t_tai.size(), 6);
      rv_in << r_in, MatX3::Zero(t_tai.size(), 3);
      MatX6 rv_out = ConvertFrameSpice(t_tai, rv_in, frame_in, frame_out);
      return rv_out.leftCols(3);
    }

    Mat6 RotOp2Mi(Real t_tai) {
      Vec6 rv_earth_icrf = GetBodyPosVel(t_tai, NaifId::SUN, NaifId::EARTH);
      Vec6 rv_moon_icrf = GetBodyPosVel(t_tai, NaifId::SUN, NaifId::MOON);

      // IAU pole
      Vec3 iau_pole = Vec3::UnitZ();
      Vec3 iau_pole_icrf = ConvertFrame(t_tai, iau_pole, Frame::MOON_PA, Frame::GCRF);

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
    ///
    /// @brief Convert the state vector from one coordinate system to another
    /// (with integer ID input)
    ///
    /// @param t_tai Epoch of the state vector
    /// @param rv_in   State vector in the original coordinate system
    /// @param frame_in Coordinate system of the original state vector
    /// @param frame_out Coordinate system of the converted state vector
    /// @return Vec6  State vector in the converted coordinate system
    ///
    Vec6 ConvertFrameSpice(Real t_tai, const Vec6& rv_in, Frame frame_in, Frame frame_out) {
      if (frame_in == frame_out) {
        return rv_in;
      }

      switch (frame_in) {
        case Frame::ITRF: {
          switch (frame_out) {
            case Frame::GCRF: {
              // Convert to GCRF
              Mat6 Rrv_itrf2gcrf = GetFrameConversionMat(t_tai, Frame::ITRF, Frame::GCRF);
              Vec6 rv_gcrf = Rrv_itrf2gcrf * rv_in;
              return rv_gcrf;
            }
            default: {  // First convert to GCRF and then to the desired frame
              Vec6 rv_gcrf = ConvertFrameSpice(t_tai, rv_in, Frame::ITRF, Frame::GCRF);
              Vec6 rv_out = ConvertFrameSpice(t_tai, rv_gcrf, Frame::GCRF, frame_out);
              return rv_out;
            }
          }
        }

        case Frame::MOON_ME: {
          switch (frame_out) {
            case Frame::MOON_PA: {  // Convert to MOON_PA
              Vec3 r_ME = rv_in.head(3);
              Vec3 v_ME = rv_in.tail(3);

              // Rotation Mat MOON_ME (in DE421) -> MOON_PA (in DE440)
              // Reference:
              // https://iopscience.iop.org/article/10.3847/1538-3881/abd414/pdf
              Mat3 B_M = RotX(-0.2785 * DEG_ARCSEC) * RotY(-78.6944 * DEG_ARCSEC)
                         * RotZ(-67.8526 * DEG_ARCSEC);
              Mat3 B_M_inv = B_M.transpose();

              Vec3 r_PA = B_M * r_ME;
              Vec3 v_PA = B_M * v_ME;
              Vec6 rv_PA;
              rv_PA << r_PA, v_PA;
              return rv_PA;
            }
            default: {  // first convert to MOON_PA and then to the desired frame
              Vec6 rv_pa = ConvertFrameSpice(t_tai, rv_in, Frame::MOON_ME, Frame::MOON_PA);
              Vec6 rv_out = ConvertFrameSpice(t_tai, rv_pa, Frame::MOON_PA, frame_out);
              return rv_out;
            }
          }
        }

        case Frame::MOON_PA: {
          switch (frame_out) {
            case Frame::MOON_ME: {
              // Rotation Mat MOON_ME (in DE421) -> MOON_PA (in DE440)
              // Reference:
              // https://iopscience.iop.org/article/10.3847/1538-3881/abd414/pdf
              Vec3 r_PA = rv_in.head(3);
              Vec3 v_PA = rv_in.tail(3);
              Mat3 B_M = RotX(-0.2785 * DEG_ARCSEC) * RotY(-78.6944 * DEG_ARCSEC)
                         * RotZ(-67.8526 * DEG_ARCSEC);
              Vec3 r_ME = B_M * r_PA;
              Vec3 v_ME = B_M * v_PA;
              Vec6 rv_me;
              rv_me << r_ME, v_ME;
              return rv_me;
            }
            case Frame::MOON_CI: {  // Convert to Moon Inertial
              Mat6 Mrot = GetFrameConversionMat(t_tai, Frame::MOON_PA, Frame::GCRF);
              Vec6 rv_mi = Mrot * rv_in;
              return rv_mi;
            }
            default: {  // first convert to MOON_CI and then to the desired frame
              Vec6 rv_mi = ConvertFrameSpice(t_tai, rv_in, Frame::MOON_PA, Frame::MOON_CI);
              Vec6 rv_out = ConvertFrameSpice(t_tai, rv_mi, Frame::MOON_CI, frame_out);
              return rv_out;
            }
          }
        }

        case Frame::GCRF: {
          switch (frame_out) {
            case Frame::ICRF: {
              Vec6 rv_icrf_ssb2e
                  = GetBodyPosVel(t_tai, NaifId::SOLAR_SYSTEM_BARYCENTER, NaifId::EARTH);
              Vec6 rv_icrf = rv_in + rv_icrf_ssb2e;
              return rv_icrf;
            }
            case Frame::ITRF: {
              Mat6 Rrv_gcrf2itrf = GetFrameConversionMat(t_tai, Frame::GCRF, Frame::ITRF);
              Vec6 rv_itrf = Rrv_gcrf2itrf * rv_in;
              return rv_itrf;
            }
            case Frame::MOON_CI: {
              Vec6 rv_icrf_m2e = GetBodyPosVel(t_tai, NaifId::MOON, NaifId::EARTH);
              Vec6 rv_mi = rv_in + rv_icrf_m2e;
              return rv_mi;
            }
            case Frame::EMR: {
              Vec6 rv_icrf_emb2e
                  = GetBodyPosVel(t_tai, NaifId::EARTH_MOON_BARYCENTER, NaifId::EARTH);
              Vec6 rv_emr = Inertial2Synodic(rv_icrf_emb2e, rv_in);
              return rv_emr;
            }
            case Frame::MOON_PA:
            case Frame::MOON_ME: {
              Vec6 rv_mi = ConvertFrameSpice(t_tai, rv_in, Frame::GCRF, Frame::MOON_CI);
              Vec6 rv_out = ConvertFrameSpice(t_tai, rv_mi, Frame::MOON_CI, frame_out);
              return rv_out;
            }
            default: assert(false && "Conversion not found");
          }
        }

        case Frame::MOON_CI: {
          switch (frame_out) {
            case Frame::GCRF: {
              Vec6 rv_icrf_e2m = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::MOON);
              Vec6 rv_gcrf = rv_in + rv_icrf_e2m;
              return rv_gcrf;
            }
            case Frame::MOON_PA: {
              Mat6 Mrot = GetFrameConversionMat(t_tai, Frame::GCRF, Frame::MOON_PA);
              Vec6 rv_pa = Mrot * rv_in;
              return rv_pa;
            }
            case Frame::MOON_ME: {
              Vec6 rv_pa = ConvertFrameSpice(t_tai, rv_in, Frame::MOON_CI, Frame::MOON_PA);
              Vec6 rv_out = ConvertFrameSpice(t_tai, rv_pa, Frame::MOON_PA, frame_out);
              return rv_out;
            }
            case Frame::MOON_OP: {
              Mat6 R_op2mi = RotOp2Mi(t_tai);
              Vec6 rv_op = R_op2mi.transpose() * rv_in;
              return rv_op;
            }
            default: {
              Vec6 rv_gcrf = ConvertFrameSpice(t_tai, rv_in, Frame::MOON_CI, Frame::GCRF);
              Vec6 rv_out = ConvertFrameSpice(t_tai, rv_gcrf, Frame::GCRF, frame_out);
              return rv_out;
            }
          }
        }

        case Frame::ICRF: {
          switch (frame_out) {
            case Frame::GCRF: {
              Vec6 rv_icrf_ssb2e
                  = GetBodyPosVel(t_tai, NaifId::SOLAR_SYSTEM_BARYCENTER, NaifId::EARTH);
              Vec6 rv_gcrf = rv_in - rv_icrf_ssb2e;
              return rv_gcrf;
            }
            default: {
              Vec6 rv_gcrf = ConvertFrameSpice(t_tai, rv_in, Frame::ICRF, Frame::GCRF);
              Vec6 rv_out = ConvertFrameSpice(t_tai, rv_gcrf, Frame::GCRF, frame_out);
              return rv_out;
            }
          }
        }

        case Frame::EMR: {
          switch (frame_out) {
            case Frame::GCRF: {
              Vec6 rv_icrf_emb2e
                  = GetBodyPosVel(t_tai, NaifId::EARTH, NaifId::EARTH_MOON_BARYCENTER);
              Vec6 rv_gcrf = Synodic2Intertial(rv_icrf_emb2e, rv_in);
              return rv_gcrf;
            }
            default: {
              Vec6 rv_gcrf = ConvertFrameSpice(t_tai, rv_in, Frame::EMR, Frame::GCRF);
              Vec6 rv_out = ConvertFrameSpice(t_tai, rv_gcrf, Frame::GCRF, frame_out);
              return rv_out;
            }
          }
        }

        case Frame::MOON_OP: {
          switch (frame_out) {
            case Frame::MOON_CI: {
              Mat6 R_op2mi = RotOp2Mi(t_tai);
              Vec6 rv_mi = R_op2mi * rv_in;
              return rv_mi;
            }
            default: {
              Vec6 rv_mi = ConvertFrameSpice(t_tai, rv_in, Frame::MOON_OP, Frame::MOON_CI);
              Vec6 rv_out = ConvertFrameSpice(t_tai, rv_mi, Frame::MOON_CI, frame_out);
              return rv_out;
            }
          }
        }
        default: {
          assert(false && "Conversion not found");
        }
      }
      assert(false && "Conversion not found");
      return Vec6::Zero();
    }

  }  // namespace spice

}  // namespace lupnt
