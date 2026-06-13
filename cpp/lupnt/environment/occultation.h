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

#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"

#pragma once

namespace lupnt {

  class Occultation {
  private:
    static constexpr double r_atmos_ = R_EARTH + 100.0e3;  // [m] Atmospheric mask
    static constexpr double r_ionos_ = R_EARTH + 965.0e3;  // [m] Ionospheric mask

  public:
    static std::map<std::string, bool> ComputeOccultation(
        Real epoch, const Vec3& r1, const Vec3& r2, Frame frame1, Frame frame2,
        const std::vector<BodyId>& bodies, const VecXd& height_atm, const VecXd& min_elevation,
        const bool use_elev_mask1, const bool use_elev_mask2);

    static std::vector<std::map<std::string, bool>> ComputeOccultation(
        Real epoch, const Mat<-1, 3>& r1, const Mat<-1, 3>& r2, Frame frame1, Frame frame2,
        const std::vector<BodyId>& bodies, const VecXd& height_atm, const VecXd& min_elevation,
        const bool use_elev_mask1, const bool use_elev_mask2);

    static std::vector<std::map<std::string, bool>> ComputeOccultation(
        const VecX& epoch, const Mat<-1, 3>& r1, const Mat<-1, 3>& r2, Frame frame1, Frame frame2,
        const std::vector<BodyId>& bodies, const VecXd& height_atm, const VecXd& min_elevation,
        const bool use_elev_mask1, const bool use_elev_mask2);
  };

}  // namespace lupnt
