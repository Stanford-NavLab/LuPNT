/**
 * @file SpiceInterface.h
 * @author Stanford NAV LAB
 * @brief  SPICE Interface functions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <string.h>

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

#include "../physics/Cheby.h"
#include "lupnt/core/Constants.h"

namespace ad = autodiff;

namespace LPT {
namespace SpiceInterface {
static segment_t *cheby_s;
static long cheby_n;
void LoadSpiceKernel(void);
void ExtractPckCoeffs(void);
Eigen::MatrixXd GetFrameConversionMatrix(ad::real et, std::string from_frame,
                                         std::string to_frame);
ad::real StringToTDB(std::string str);
ad::real StringToTAI(std::string str);
std::string TAItoStringUTC(ad::real tai, int prec);
std::string TDBtoStringUTC(ad::real tdb, int prec);
ad::real ConvertTime(ad::real t, std::string from_time_type,
                     std::string to_time_type);
ad::VectorXreal GetBodyPosVel(const ad::real tai, int center, int target);
Eigen::Vector3d GetBodyPos(std::string targetName, ad::real epoch,
                           std::string refFrame, std::string obsName,
                           std::string abCorrection);

}  // namespace  SpiceInterface

}  // namespace LPT
