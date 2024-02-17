/**
 * @file CoordConverter.cpp
 * @author Stanford NAV LAB
 * @brief Coordinate conversion functions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "coord_converter.h"

#include <filesystem>

#include "cheby.h"
#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/orbit_state.h"
#include "lupnt/physics/orbit_state_utils.h"
#include "spice_interface.h"

using namespace lupnt;
using namespace SpiceInterface;

CartesianOrbitState CoordConverter::Convert(real epoch,
                                            const CartesianOrbitState &state_in,
                                            CoordSystem coord_sys_out) {
  Vector6 rv_in = state_in.GetVector();
  Vector6 rv_out =
      Convert(epoch, rv_in, state_in.GetCoordSystem(), coord_sys_out);
  return CartesianOrbitState(rv_out, coord_sys_out);
}

Vector3 CoordConverter::Convert(real epoch, Vector3 rv_in,
                                CoordSystem coord_sys_in,
                                CoordSystem coord_sys_out) {
  Vector6 rv_in_6;
  rv_in_6 << rv_in, Vector3::Zero();
  Vector6 rv_out_6 = Convert(epoch, rv_in_6, coord_sys_in, coord_sys_out);
  return rv_out_6.head(3);
}

/**
 * @brief Convert the state vector from one coordinate system to another (with
 * integer ID input)
 *
 * @param epoch Epoch of the state vector
 * @param rv_in   State vector in the original coordinate system
 * @param coord_sys_in Coordinate system of the original state vector
 * @param coord_sys_out Coordinate system of the converted state vector
 * @return Vector6  State vector in the converted coordinate system
 */
Vector6 CoordConverter::Convert(real epoch, Vector6 rv_in,
                                CoordSystem coord_sys_in,
                                CoordSystem coord_sys_out) {
  if (coord_sys_in == coord_sys_out) {
    return rv_in;
  }

  switch (coord_sys_in) {
    case ITRF: {
      switch (coord_sys_out) {
        case GCRF: {
          // Convert to GCRF
          Matrix6 Rrv_gcrf_itrf = ComputeITRFtoGCRF(epoch);
          Vector6 rv_gcrf = Rrv_gcrf_itrf * rv_in;
          return rv_gcrf;
        }
        default: {  // First convert to GCRF and then to the desired frame
          Vector6 rv_gcrf = Convert(epoch, rv_in, ITRF, GCRF);
          Vector6 rv_out = Convert(epoch, rv_gcrf, GCRF, coord_sys_out);
          return rv_out;
        }
      }
    }

    case ME: {
      switch (coord_sys_out) {
        case PA: {  // Convert to PA
          Vector3 r_ME = rv_in.head(3);
          Vector3 v_ME = rv_in.tail(3);

          // Rotation Matrix ME (in DE421) -> PA (in DE440)
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
        default: {  // first convert to PA and then to the desired frame
          Vector6 rv_pa = Convert(epoch, rv_in, ME, PA);
          Vector6 rv_out = Convert(epoch, rv_pa, PA, coord_sys_out);
          return rv_out;
        }
      }
    }

    case PA: {
      switch (coord_sys_out) {
        case ME: {
          // Rotation Matrix ME (in DE421) -> PA (in DE440)
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
        case MI: {  // Convert to Moon Inertial
          Matrix6 Mrot = GetFrameConversionMatrix(epoch, "MOON_PA", "J2000");
          Vector6 rv_mi = Mrot * rv_in;
          return rv_mi;
        }
        default: {  // first convert to MI and then to the desired frame
          Vector6 rv_mi = Convert(epoch, rv_in, PA, MI);
          Vector6 rv_out = Convert(epoch, rv_mi, MI, coord_sys_out);
          return rv_out;
        }
      }
    }

    case GCRF: {
      switch (coord_sys_out) {
        case ICRF: {
          Vector6 rv_icrf_ssb_earth = GetBodyPosVel(
              epoch, NaifId::SOLAR_SYSTEM_BARYCENTER, NaifId::EARTH);
          Vector6 rv_icrf = rv_in + rv_icrf_ssb_earth;
          return rv_icrf;
        }
        case ITRF: {
          Matrix6 Rrv_itrf_gcrf = ComputeGCRFtoITRF(epoch);
          Vector6 rv_itrf = Rrv_itrf_gcrf * rv_in;
          return rv_itrf;
        }
        case MI: {
          Vector6 rv_icrf_moon_earth =
              GetBodyPosVel(epoch, NaifId::MOON, NaifId::EARTH);
          Vector6 rv_mi = rv_in + rv_icrf_moon_earth;
          return rv_mi;
        }
        case EMR: {
          Vector6 rv_icrf_emb_earth = GetBodyPosVel(
              epoch, NaifId::EARTH_MOON_BARYCENTER, NaifId::EARTH);
          Vector6 rv_emr = InertialToRtn(rv_icrf_emb_earth, rv_in);
          return rv_emr;
        }
        case PA:
        case ME: {
          Vector6 rv_mi = Convert(epoch, rv_in, GCRF, MI);
          Vector6 rv_out = Convert(epoch, rv_mi, MI, coord_sys_out);
          return rv_out;
        }
        default:
          assert(false && "Conversion not found");
      }
    }

    case MI: {
      switch (coord_sys_out) {
        case GCRF: {
          Vector6 rv_icrf_earth_moon =
              GetBodyPosVel(epoch, NaifId::EARTH, NaifId::MOON);
          Vector6 rv_gcrf = rv_in + rv_icrf_earth_moon;
          return rv_gcrf;
        }
        case PA: {
          Matrix6 Mrot = GetFrameConversionMatrix(epoch, "J2000", "MOON_PA");
          Vector6 rv_pa = Mrot * rv_in;
          return rv_pa;
        }
        case ME: {
          Vector6 rv_pa = Convert(epoch, rv_in, MI, PA);
          Vector6 rv_out = Convert(epoch, rv_pa, PA, coord_sys_out);
          return rv_out;
        }
        case OP: {
          Matrix6 R_op2mi = ComputeOpToMi(epoch);
          Vector6 rv_op = R_op2mi.transpose() * rv_in;
          return rv_op;
        }
        default: {
          Vector6 rv_gcrf = Convert(epoch, rv_in, MI, GCRF);
          Vector6 rv_out = Convert(epoch, rv_gcrf, GCRF, coord_sys_out);
          return rv_out;
        }
      }
    }

    case ICRF: {
      switch (coord_sys_out) {
        case GCRF: {
          Vector6 rv_icrf_ssb_earth = GetBodyPosVel(
              epoch, NaifId::SOLAR_SYSTEM_BARYCENTER, NaifId::EARTH);
          Vector6 rv_gcrf = rv_in - rv_icrf_ssb_earth;
          return rv_gcrf;
        }
        default: {
          Vector6 rv_gcrf = Convert(epoch, rv_in, ICRF, GCRF);
          Vector6 rv_out = Convert(epoch, rv_gcrf, GCRF, coord_sys_out);
          return rv_out;
        }
      }
    }

    case EMR: {
      switch (coord_sys_out) {
        case GCRF: {
          Vector6 rv_icrf_emb_earth = GetBodyPosVel(
              epoch, NaifId::EARTH, NaifId::EARTH_MOON_BARYCENTER);
          Vector6 rv_gcrf = RtnToInertial(rv_icrf_emb_earth, rv_in);
          return rv_gcrf;
        }
        default: {
          Vector6 rv_gcrf = Convert(epoch, rv_in, EMR, GCRF);
          Vector6 rv_out = Convert(epoch, rv_gcrf, GCRF, coord_sys_out);
          return rv_out;
        }
      }
    }

    case OP: {
      switch (coord_sys_out) {
        case MI: {
          Matrix6 R_op2mi = ComputeOpToMi(epoch);
          Vector6 rv_mi = R_op2mi * rv_in;
          return rv_mi;
        }
        default: {
          Vector6 rv_mi = Convert(epoch, rv_in, OP, MI);
          Vector6 rv_out = Convert(epoch, rv_mi, MI, coord_sys_out);
          return rv_out;
        }
      }
    }

    default: {
      assert(false && "Conversion not found");
    }
  }
}

/**
 * @brief Compute the rotation matrix that captures the ‘precession-nutation’
 * component of the transformation between body-fixed and inertial coordinate
 * systems
 *
 * @param epoch  tai in JD
 * @return Matrix3 3x3 rotation matrix
 * @ref NASA/TP–20220014814 page 31-36, GMAT Math Spec page 19, GMAT
 * ITRFAxes::CalculateRotationMatrix
 */
Matrix6 CoordConverter::ComputeITRFtoGCRF(real tai) {
  // Get the rotation matrix using SPICE
  auto tdb_s = ConvertTime(tai, "TAI", "TDB");
  double et = tdb_s.val();  // this cuts of the relationship between t and
                            // Mrot temporaly
  Matrix6 Mrot = GetFrameConversionMatrix(et, "ITRF93", "J2000");
  return Mrot;
}

Matrix6 CoordConverter::ComputeGCRFtoITRF(real tai) {
  // Get the rotation matrix using SPICE
  auto tdb_s = ConvertTime(tai, "TAI", "TDB");
  double et = tdb_s.val();  // this cuts of the relationship between t and
                            // Mrot temporaly
  Matrix6 Mrot = GetFrameConversionMatrix(et, "J2000", "ITRF93");
  return Mrot;
}

Matrix6 CoordConverter::ComputeOpToMi(real epoch) {
  Vector6 rv_earth_icrf = GetBodyPosVel(epoch, NaifId::SUN, NaifId::EARTH);
  Vector6 rv_moon_icrf = GetBodyPosVel(epoch, NaifId::SUN, NaifId::MOON);

  // Moon axis
  MatrixX iau_moon2icrf;
  iau_moon2icrf =
      GetFrameConversionMatrix(epoch, "IAU_MOON", "J2000").block(0, 0, 3, 3);

  // IAU pole
  Vector3 iau_pole{0, 1, 0};
  Vector3 iau_pole_icrf = iau_moon2icrf * iau_pole;

  // OP unit vectors
  Vector3 dr = rv_earth_icrf.head(3) - rv_moon_icrf.head(3);
  Vector3 dv = rv_earth_icrf.tail(3) - rv_moon_icrf.tail(3);
  Vector3 z_op = dr.cross(dv).normalized();
  Vector3 x_op = iau_pole_icrf.cross(z_op).normalized();
  Vector3 y_op = z_op.cross(x_op).normalized();

  // create rotation matrix from OP to MI
  Matrix3 R_op2mi;
  R_op2mi << x_op, y_op, z_op;

  Matrix6 R_op2mi_tot = Matrix6::Zero();
  R_op2mi_tot.topLeftCorner(3, 3) = R_op2mi;
  R_op2mi_tot.bottomRightCorner(3, 3) = R_op2mi;

  return R_op2mi_tot.cast<double>();
}
