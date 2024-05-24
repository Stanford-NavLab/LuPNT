/**
 * @file occultation.h
 * @author Stanford NAV LAB
 * @brief Signal blockage model
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <cmath>
#include <map>

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/frame_converter.h"

#pragma once

namespace lupnt {

class Occultation {
 private:
  static constexpr double min_elevation_ = 10.0 * RAD_PER_DEG;
  static constexpr double r_atmos_ = R_EARTH + 100.0;  // atmospheric mask
  static constexpr double r_ionos_ = R_EARTH + 965.0;  // ionospheric mask

 public:
  /**
   * @brief Compute occultation between a tx and a rx
   *
   * @param tx_eci   The tx position in the Earth-centered inertial
   * @param tx_mci   The tx position in the Moon-centered inertial
   * @param rx_eci      The rx position in the Earth-centered inertial frame
   * @param rx_mci      The rx position in the Moon-centered inertial frame
   * @param tx_planet    The planet surface the transmitter is on
   * @return std::map<std::string, bool>
   */
  static std::map<std::string, bool> ComputeOccultation(
      const Vector3d tx_eci, const Vector3d tx_mci, const Vector3d rx_eci,
      const Vector3d rx_mci, const std::string tx_planet);

  // VectorX = func(real, Vector3, Vector3, ...)
  static VectorX ComputeOccultation(real epoch, const Vector3& r1,
                                    const Vector3& r2, Frame cs1, Frame cs2,
                                    const std::vector<NaifId>& bodies);

  // MatrixX = func(real, Matrix<-1, 3>, Matrix<-1, 3>, ...)
  static MatrixX ComputeOccultation(real epoch, const Matrix<-1, 3>& r1,
                                    const Matrix<-1, 3>& r2, Frame cs1,
                                    Frame cs2,
                                    const std::vector<NaifId>& bodies);

  // MatrixX = func(VectorX, Matrix<-1, 3>, Matrix<-1, 3>, ...)
  static MatrixX ComputeOccultation(const VectorX& epoch,
                                    const Matrix<-1, 3>& r1,
                                    const Matrix<-1, 3>& r2, Frame cs1,
                                    Frame cs2,
                                    const std::vector<NaifId>& bodies);
};
}  // namespace lupnt
