/**
 * @file spice_interface.h
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

#include "lupnt/core/constants.h"
#include "lupnt/physics/cheby.h"

namespace lupnt {
namespace SpiceInterface {

static segment_t *cheby_s;
static long cheby_n;
void LoadSpiceKernel(void);
void ExtractPckCoeffs(void);
Matrix6d GetFrameConversionMatrix(real et, std::string from_frame,
                                  std::string to_frame);
real StringToTDB(std::string str);
real StringToTAI(std::string str);
std::string TAItoStringUTC(real tai, int prec);
std::string TDBtoStringUTC(real tdb, int prec);
real ConvertTime(real t, std::string from_time_type, std::string to_time_type);
Vector6 GetBodyPosVel(const real tai, NaifId center, NaifId target);
Matrix<-1, 6> GetBodyPosVel(const VectorX &tai, NaifId center, NaifId target);
Vector3d GetBodyPos(std::string targetName, real epoch, std::string refFrame,
                    std::string obsName, std::string abCorrection);

}  // namespace SpiceInterface

}  // namespace lupnt
