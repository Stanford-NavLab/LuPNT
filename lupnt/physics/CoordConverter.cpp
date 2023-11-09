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

#include "CoordConverter.h"

#include <autodiff/forward/real/eigen.hpp>
#include <filesystem>

#include "../core/Constants.h"
#include "../numerics/MathUtils.h"
#include "../physics/SpiceInterface.h"
#include "Cheby.h"

using namespace LPT;
using namespace SpiceInterface;

const std::string CoordConverter::COORD_SYSTEM_TEXT[CoordSystemCount] = {
    "ITRF", "GCRF", "ICRF", "SER", "GSE", "EME", "MOD",
    "TOD",  "EMR",  "MI",   "PA",  "ME",  "RTN"};

CoordSystem CoordConverter::GetCoordTypeID(const std::string &str) {
  for (int i = 0; i < CoordSystemCount; i++) {
    if (str == COORD_SYSTEM_TEXT[i]) return static_cast<CoordSystem>(i);
  }
  throw std::runtime_error(
      "CoordConverter::GetCoordTypeID: Unknown coordinate system ID");
}

/**
 * @brief Convert the state vector from one coordinate system to another (with
 * string ID input)
 *
 * @param rv_in  State vector
 * @param epoch  TAI epoch [s]
 * @param coord_sys_in   Coordinate system ID
 * @param coord_sys_out     Coordinate system ID
 * @return ad::VectorXreal      State vector
 */
ad::VectorXreal CoordConverter::Convert(const ad::VectorXreal rv_in,
                                        const ad::real epoch,
                                        const std::string coord_sys_in,
                                        const std::string coord_sys_out) {
  return Convert(rv_in, epoch, GetCoordTypeID(coord_sys_in),
                 GetCoordTypeID(coord_sys_out));
}

/**
 * @brief Convert the state vector from one coordinate system to another (with
 * integer ID input)
 *
 * @param rv_in   State vector in the original coordinate system
 * @param epoch Epoch of the state vector
 * @param coord_sys_in Coordinate system of the original state vector
 * @param coord_sys_out Coordinate system of the converted state vector
 * @return ad::VectorXreal  State vector in the converted coordinate system
 */
ad::VectorXreal CoordConverter::Convert(const ad::VectorXreal rv_in,
                                        const ad::real epoch,
                                        const CoordSystem coord_sys_in,
                                        const CoordSystem coord_sys_out) {
  if (coord_sys_in == coord_sys_out) return rv_in;

  ad::VectorXreal rv_out(6);

  if (coord_sys_in == ITRF)  /// Earth fixed frame
  {
    if (coord_sys_out == GCRF) {  // Earth fixed -> Inertial Frame
      // Convert to GCRF
      ad::MatrixXreal Rrv_GCRF_ITRF(6, 6);
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
      ad::Vector3real r_ME = rv_in.head(3);
      ad::Vector3real v_ME = rv_in.tail(3);

      // Rotation Matrix ME (in DE421) -> PA (in DE440)
      // Reference:
      // https://iopscience.iop.org/article/10.3847/1538-3881/abd414/pdf
      ad::MatrixXreal B_M = R1(-0.2785 * DEG_PER_ARCSEC) *
                            R2(-78.6944 * DEG_PER_ARCSEC) *
                            R3(-67.8526 * DEG_PER_ARCSEC);
      ad::MatrixXreal B_M_inv = B_M.transpose();

      ad::Vector3real r_PA = B_M * r_ME;
      ad::Vector3real v_PA = B_M * v_ME;
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
      ad::Vector3real r_PA = rv_in.head(3);
      ad::Vector3real v_PA = rv_in.tail(3);
      ad::MatrixXreal B_M = R1(-0.2785 * DEG_PER_ARCSEC) *
                            R2(-78.6944 * DEG_PER_ARCSEC) *
                            R3(-67.8526 * DEG_PER_ARCSEC);
      ad::Vector3real r_ME = B_M * r_PA;
      ad::Vector3real v_ME = B_M * v_PA;
      rv_out << r_ME, v_ME;
      return rv_out;
    } else if (coord_sys_out == MI)  // Convert to Moon Inertial
    {
      ad::MatrixXreal Mrot =
          GetFrameConversionMatrix(epoch, "MOON_PA", "J2000");
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
    ad::Vector3real r_GCRF_E_p = rv_in.head(3);
    ad::Vector3real v_GCRF_E_p = rv_in.tail(3);

    switch (coord_sys_out) {
      case ICRF: {
        // Get position and velocity of the Solar System Barycenter (SSB) with
        // respect to Earth (E) in ICRF frame
        ad::VectorXreal rv_ICRF_SSB_E = GetBodyPosVel(epoch, 399, 0);
        rv_out = rv_in + rv_ICRF_SSB_E;
        return rv_out;
      }
      case ITRF: {
        // Convert to GCRF
        ad::MatrixXreal Rrv_GCRF_ITRF(6, 6);
        Rrv_GCRF_ITRF = ComputeITRFtoGCRF(epoch);
        rv_out = Rrv_GCRF_ITRF.transpose() * rv_in;

        return rv_out;
      }
      case MI: {
        // Get position and velocity of the SMoon with respect to Earth (E) in
        // ICRF frame
        ad::VectorXreal rv_ICRF_Moon_E = GetBodyPosVel(epoch, 399, 301);
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
      ad::VectorXreal rv_ICRF_E_Moon = GetBodyPosVel(epoch, 301, 399);
      rv_out = rv_in + rv_ICRF_E_Moon;  // (sc - M) + (Moon - E) = (sc - Moon)
      return rv_out;
    } else if (coord_sys_out == PA) {  //  Convert to PA
      // convert TAI to TDB past J2000
      ad::MatrixXreal Mrot =
          GetFrameConversionMatrix(epoch, "J2000", "MOON_PA");
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
 * @return ad::Matrix3real 3x3 rotation matrix
 * @ref NASA/TP–20220014814 page 31-36, GMAT Math Spec page 19, GMAT
 * ITRFAxes::CalculateRotationMatrix
 */
ad::MatrixXreal CoordConverter::ComputeITRFtoGCRF(const ad::real tai) {
  // Get the rotation matrix using SPICE
  ad::real tdb_s = ConvertTime(tai, "TAI", "TDB");
  double et =
      tdb_s
          .val();  // this cuts of the relationship between t and Mrot temporaly
  ad::MatrixXreal Mrot(6, 6);
  Mrot = GetFrameConversionMatrix(et, "ITRF93", "J2000");

  return Mrot;
}
