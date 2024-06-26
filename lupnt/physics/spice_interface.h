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
Mat6d GetFrameConversionMat(real t_tai, Frame from_frame, Frame to_frame);

real DateToModifiedJulianDate(int year, int month, int day, int hour,
                              int minute, real second);
std::tuple<int, int, int, int, int, real> ModifiedJulianDateToDate(real mjd);

real StringToTDB(std::string str);
real StringToTAI(std::string str);

std::string TAItoStringUTC(real t_tai, int prec);
std::string TDBtoStringUTC(real t_tdb, int prec);

real EarthRotationAngle(real mjd_ut1);

real GreenwichMeanSiderealTime(real mjd_ut1);
real UTCtoUT1(real mjd_utc);
real ConvertTime(real t, std::string from_time, std::string to_time);

Vec6 GetBodyPosVel(const real t_tai, NaifId center, NaifId target, Frame frame);
Mat<-1, 6> GetBodyPosVel(const VecX &t_tai, NaifId center, NaifId target,
                         Frame frame);
Vec3d GetBodyPos(NaifId target, real t_tai, Frame refFrame, NaifId obs,
                 std::string abCorrection);

}  // namespace lupnt
