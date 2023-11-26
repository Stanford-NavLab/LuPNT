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

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
#include "cheby.h"
#include "spice_interface.h"

using namespace lupnt;
using namespace SpiceInterface;

/**
 * @brief Convert the state vector from one coordinate system to another (with
 * integer ID input)
 *
 * @param rv_in   State vector in the original coordinate system
 * @param epoch Epoch of the state vector
 * @param coord_sys_in Coordinate system of the original state vector
 * @param coord_sys_out Coordinate system of the converted state vector
 * @return Vector6  State vector in the converted coordinate system
 */
Vector6 CoordConverter::Convert(Vector6 rv_in, real epoch,
                                    CoordSystem coord_sys_in,
                                    CoordSystem coord_sys_out) {
  if (coord_sys_in == coord_sys_out) return rv_in;

  Vector6 rv_out;

  if (coord_sys_in == ITRF)  /// Earth fixed frame
  {
    if (coord_sys_out == GCRF) {  // Earth fixed -> Inertial Frame
      // Convert to GCRF
      Matrix6 Rrv_GCRF_ITRF;
      Rrv_GCRF_ITRF = ComputeITRFtoGCRF(epoch);
      rv_out = Rrv_GCRF_ITRF * rv_in;
      return rv_out;
    } else  // First convert to GCRF and then to the desired frame
    {
      rv_out = Convert(rv_in, epoch, ITRF, GCRF);
      rv_out = Convert(rv_out, epoch, GCRF, coord_sys_out);
      return rv_out;
    }
  }

  if (coord_sys_in == ME) {
    if (coord_sys_out == PA) {
      // Convert to PA
      Vector3 r_ME = rv_in.head(3);
      Vector3 v_ME = rv_in.tail(3);

      // Rotation Matrix ME (in DE421) -> PA (in DE440)
      // Reference:
      // https://iopscience.iop.org/article/10.3847/1538-3881/abd414/pdf
      Matrix3 B_M = R1(-0.2785 * DEG_PER_ARCSEC) *
                        R2(-78.6944 * DEG_PER_ARCSEC) *
                        R3(-67.8526 * DEG_PER_ARCSEC);
      Matrix3 B_M_inv = B_M.transpose();

      Vector3 r_PA = B_M * r_ME;
      Vector3 v_PA = B_M * v_ME;
      rv_out << r_PA, v_PA;
      return rv_out;
    } else  // first convert to PA and then to the desired frame
    {
      rv_out = Convert(rv_in, epoch, ME, PA);
      rv_out = Convert(rv_out, epoch, PA, coord_sys_out);
      return rv_out;
    }
  }

  if (coord_sys_in == PA) {
    if (coord_sys_out == ME) {
      // Rotation Matrix ME (in DE421) -> PA (in DE440)
      // Reference:
      // https://iopscience.iop.org/article/10.3847/1538-3881/abd414/pdf
      Vector3 r_PA = rv_in.head(3);
      Vector3 v_PA = rv_in.tail(3);
      Matrix3 B_M = R1(-0.2785 * DEG_PER_ARCSEC) *
                        R2(-78.6944 * DEG_PER_ARCSEC) *
                        R3(-67.8526 * DEG_PER_ARCSEC);
      Vector3 r_ME = B_M * r_PA;
      Vector3 v_ME = B_M * v_PA;
      rv_out << r_ME, v_ME;
      return rv_out;
    } else if (coord_sys_out == MI)  // Convert to Moon Inertial
    {
      Matrix6 Mrot = GetFrameConversionMatrix(epoch, "MOON_PA", "J2000");
      rv_out = Mrot * rv_in;
      return rv_out;
    } else {  // first convert to MI and then to the desired frame
      rv_out = Convert(rv_in, epoch, PA, MI);
      rv_out = Convert(rv_out, epoch, MI, coord_sys_out);
      return rv_out;
    }
  }

  if (coord_sys_in == GCRF)  // Earth centered Inertial Frame
  {
    Vector3 r_GCRF_E_p = rv_in.head(3);
    Vector3 v_GCRF_E_p = rv_in.tail(3);

    switch (coord_sys_out) {
      case ICRF: {
        // Get position and velocity of the Solar System Barycenter (SSB) with
        // respect to Earth (E) in ICRF frame
        VectorX rv_ICRF_SSB_E = GetBodyPosVel(epoch, 399, 0);
        rv_out = rv_in + rv_ICRF_SSB_E;
        return rv_out;
      }
      case ITRF: {
        // Convert to GCRF
        Matrix6 Rrv_GCRF_ITRF;
        Rrv_GCRF_ITRF = ComputeITRFtoGCRF(epoch);
        rv_out = Rrv_GCRF_ITRF.transpose() * rv_in;

        return rv_out;
      }
      case MI: {
        // Get position and velocity of the SMoon with respect to Earth (E) in
        // ICRF frame
        Vector6 rv_ICRF_Moon_E = GetBodyPosVel(epoch, 399, 301);
        rv_out = rv_in + rv_ICRF_Moon_E;  // (sc - E) + (E - Moon) = (sc - Moon)
        return rv_out;
      }
      default:
        throw std::runtime_error("Not implemented yet");
    }
  }

  if (coord_sys_in == MI) {
    if (coord_sys_out == GCRF) {
      // Get position and velocity of the SMoon with respect to Earth (E) in
      // ICRF frame
      Vector6 rv_ICRF_E_Moon = GetBodyPosVel(epoch, 301, 399);
      rv_out = rv_in + rv_ICRF_E_Moon;  // (sc - M) + (Moon - E) = (sc - Moon)
      return rv_out;
    } else if (coord_sys_out == PA) {  //  Convert to PA
      // convert TAI to TDB past J2000
      Matrix6 Mrot = GetFrameConversionMatrix(epoch, "J2000", "MOON_PA");
      rv_out = Mrot * rv_in;
      return rv_out;
    } else if (coord_sys_out == ME) {  // first convert to PA and then to ME
      rv_out = Convert(rv_in, epoch, MI, PA);
      rv_out = Convert(rv_out, epoch, PA, coord_sys_out);
      return rv_out;
    } else {  // first convert to ME and then to the desired frame
      rv_out = Convert(rv_in, epoch, MI, GCRF);
      rv_out = Convert(rv_out, epoch, GCRF, coord_sys_out);
      return rv_out;
    }
  }

  throw std::runtime_error("Conversion not implemented");
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
  double et =
      tdb_s
          .val();  // this cuts of the relationship between t and Mrot temporaly
  Matrix6 Mrot = GetFrameConversionMatrix(et, "ITRF93", "J2000");

  return Mrot;
}

/**
 * @brief Compute the rotation matrix corresponding to a simple rotation about
 * the first axis in a system
 *
 * @param phi angle of rotation
 * @return Matrix3 3x3 rotation matrix
 * @ref NASA/TP–20220014814 page 30
 */
Matrix3 CoordConverter::R1(real phi) {
  real c = cos(phi);
  real s = sin(phi);
  Matrix3 R1{
      {1.0, 0.0, 0.0},
      {0.0, c, s},
      {0.0, -s, c},
  };
  return R1;
}

/**
 * @brief Compute the rotation matrix corresponding to a simple rotation about
 * the second axis in a system
 *
 * @param phi angle of rotation
 * @return Matrix3 3x3 rotation matrix
 * @ref NASA/TP–20220014814 page 30
 */
Matrix3 CoordConverter::R2(real phi) {
  real c = cos(phi);
  real s = sin(phi);
  Matrix3 R2{
      {c, 0.0, -s},
      {0.0, 1.0, 0.0},
      {s, 0.0, c},
  };
  return R2;
}

/**
 * @brief Compute the rotation matrix corresponding to a simple rotation about
 * the third axis in a system
 *
 * @param phi angle of rotation
 * @return Matrix3 3x3 rotation matrix
 * @ref NASA/TP–20220014814 page 30
 */
Matrix3 CoordConverter::R3(real phi) {
  real c = cos(phi);
  real s = sin(phi);
  Matrix3 R3{
      {c, s, 0.0},
      {-s, c, 0.0},
      {0.0, 0.0, 1.0},
  };
  return R3;
}

/**
 * @brief Compute the rotation matrix specified by the input skew vector.
 *
 * @param v
 * @return Matrix3
 */
Matrix3 CoordConverter::Skew(Vector3 x)
{
  Matrix3 skew{
      {0.0, -x(2), x(1)},
      {x(2), 0.0, -x(0)},
      {-x(1), x(0), 0.0},
  };
  return skew;
}