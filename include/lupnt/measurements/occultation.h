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
    static constexpr double min_elevation_ = 10.0 * RAD;
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
    static std::map<std::string, bool> ComputeOccultation(const Vec3d tx_eci, const Vec3d tx_mci,
                                                          const Vec3d rx_eci, const Vec3d rx_mci,
                                                          const std::string tx_planet);

    // VecX = func(real, Vec3, Vec3, ...)
    static VecX ComputeOccultation(Real epoch, const Vec3& r1, const Vec3& r2, Frame cs1, Frame cs2,
                                   const std::vector<NaifId>& bodies);

    // MatX = func(real, Mat<-1, 3>, Mat<-1, 3>, ...)
    static MatX ComputeOccultation(Real epoch, const Mat<-1, 3>& r1, const Mat<-1, 3>& r2,
                                   Frame cs1, Frame cs2, const std::vector<NaifId>& bodies);

    // MatX = func(VecX, Mat<-1, 3>, Mat<-1, 3>, ...)
    static MatX ComputeOccultation(const VecX& epoch, const Mat<-1, 3>& r1, const Mat<-1, 3>& r2,
                                   Frame cs1, Frame cs2, const std::vector<NaifId>& bodies);
  };
}  // namespace lupnt
