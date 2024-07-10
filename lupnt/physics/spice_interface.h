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
#include <tuple>

#include "lupnt/core/constants.h"
#include "lupnt/physics/cheby.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {

// map from CoordSytem to string
const std::map<Frame, std::string> frametem_string = {
    {Frame::ITRF, "ITRF93"},
    {Frame::GCRF, "J2000"},
    {Frame::MOON_PA, "MOON_PA"},
    {Frame::MOON_CI, "J2000"},
};

static segment_t *cheby_s;
static long cheby_n;
void LoadSpiceKernel(void);
void ExtractPckCoeffs(void);
Mat6d GetFrameConversionMat(Real t_tai, Frame from_frame, Frame to_frame);

Real String2TDB(std::string str);
Real String2TAI(std::string str);

std::string TAItoStringUTC(Real t_tai, int prec);
std::string TDBtoStringUTC(Real t_tdb, int prec);

Real ConvertTime(Real t, std::string from_time, std::string to_time);

Vec6 GetBodyPosVel(const Real t_tai, NaifId center, NaifId target);
Mat<-1, 6> GetBodyPosVel(const VecX &t_tai, NaifId center, NaifId target);
Vec3d GetBodyPosSpice(NaifId target, Real t_tai, Frame refFrame, NaifId obs,
                      std::string abCorrection);

}  // namespace lupnt
