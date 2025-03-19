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
     * @brief Compute occultation between a tx and a rx for Lunar satellites
     *
     * @param tx_eci   The tx position in the Earth-centered inertial
     * @param tx_mci   The tx position in the Moon-centered inertial
     * @param rx_eci      The rx position in the Earth-centered inertial frame
     * @param rx_mci      The rx position in the Moon-centered inertial frame
     * @param tx_planet    The planet surface the transmitter is on
     * @return std::map<std::string, bool>
     */
    static std::map<std::string, bool> ComputeOccultationGnss(const Vec3d tx_eci,
                                                              const Vec3d tx_mci,
                                                              const Vec3d rx_eci,
                                                              const Vec3d rx_mci,
                                                              const std::string tx_planet);

    /**
     * @brief Compute occultation between a tx and a rx for a list of planets
     *
     * @param epoch  The epoch in TAI
     * @param r1   The tx position
     * @param r2   The tx position
     * @param cs1  The frame of the tx
     * @param cs2  The frame of the rx
     * @param bodies  The list of planets to check for occultation (naif ids)
     * @param atm_h  The atmospheric heights of the planets
     * @return std::map<string, bool>  A map of the planets and their occultation,
     * ["all"] is the total occultation
     */
    static std::map<std::string, bool> ComputeOccultation(Real epoch, const Vec3& r1,
                                                          const Vec3& r2, Frame cs1, Frame cs2,
                                                          const std::vector<NaifId>& bodies,
                                                          const VecXd& atm_h);

    /**
     * @brief  Compute occultation between a tx and a rx for a list of planets
     * (vectorized, epoch fixed)
     *
     * @param epoch  The epoch in TAI
     * @param r1   The tx position
     * @param r2   The rx position
     * @param cs1  The frame of the tx
     * @param cs2  The frame of the rx
     * @param bodies   The list of planets to check for occultation (naif ids)
     * @return std::vector<std::map<string, bool>>
     */
    static std::vector<std::map<std::string, bool>> ComputeOccultation(
        Real epoch, const Mat<-1, 3>& r1, const Mat<-1, 3>& r2, Frame cs1, Frame cs2,
        const std::vector<NaifId>& bodies, const VecXd& atm_h);

    /**
     * @brief  Compute occultation between a tx and a rx for a list of planets
     * (vectorized, epoch vector)
     *
     * @param epoch  The epoch vector in TAI
     * @param r1   The tx position
     * @param r2   The rx position
     * @param cs1  The frame of the tx
     * @param cs2  The frame of the rx
     * @param bodies   The list of planets to check for occultation (naif ids)
     * @return std::vector<std::map<string, bool>>
     */
    static std::vector<std::map<std::string, bool>> ComputeOccultation(
        const VecX& epoch, const Mat<-1, 3>& r1, const Mat<-1, 3>& r2, Frame cs1, Frame cs2,
        const std::vector<NaifId>& bodies, const VecXd& atm_h);
  };
}  // namespace lupnt
