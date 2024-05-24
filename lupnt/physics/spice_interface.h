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

#include <map>

#include "lupnt/core/constants.h"
#include "lupnt/physics/cheby.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {
namespace SpiceInterface {

// map from CoordSytem to string
const std::map<Frame, std::string> frametem_string = {
    {Frame::ITRF, "ITRF93"},
    {Frame::GCRF, "J2000"},
    {Frame::PA, "MOON_PA"},
};

static segment_t *cheby_s;
static long cheby_n;
void LoadSpiceKernel(void);
void ExtractPckCoeffs(void);
Matrix6d GetFrameConversionMatrix(real t_tai, Frame from_frame, Frame to_frame);

real StringToTDB(std::string str);
real StringToTAI(std::string str);

std::string TAItoStringUTC(real t_tai, int prec);
std::string TDBtoStringUTC(real t_tdb, int prec);

real ConvertTime(real t, std::string from_time, std::string to_time);

Vector6 GetBodyPosVel(const real t_tai, NaifId center, NaifId target);
Matrix<-1, 6> GetBodyPosVel(const VectorX &t_tai, NaifId center, NaifId target);
Vector3d GetBodyPos(NaifId target, real t_tai, Frame refFrame, NaifId obs,
                    std::string abCorrection);

}  // namespace SpiceInterface

}  // namespace lupnt
