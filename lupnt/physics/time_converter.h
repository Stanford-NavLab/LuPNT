#pragma once

#include <string>
#include <tuple>

#include "lupnt/core/constants.h"
#include "lupnt/numerics/graphs.h"
#include "lupnt/numerics/vector_macros.h"

namespace lupnt {

extern std::map<std::pair<std::string, std::string>, std::function<Real(Real)>>
    time_conversions;

Real ConvertT(Real t, const std::string& from, const std::string& to);
VecX ConvertT(VecX t, const std::string& from, const std::string& to);

Real UTCtoUT1(Real t_utc);
Real UT1toUTC(Real t_ut1);

Real TAItoUTC(Real t_tai);
Real UTCtoTAI(Real t_utc);

Real TAItoTT(Real t_tai);
Real TTtoTAI(Real t_tt);

Real TCGtoTT(Real t_tcg);
Real TTtoTCG(Real t_tt);

Real TTtoTDB(Real t_tt);
Real TDBtoTT(Real t_tdb);

Real TAItoGPS(Real t_tai);
Real GPStoTAI(Real t_gps);

Real TCBtoTDB(Real t_tcb);
Real TTtoTCB(Real t_tdb);

Real MJDtoTime(Real mjd);
Real TimeToMJD(Real t);

Real JDtoTime(Real jd);
Real TimeToJD(Real t);

Real EarthRotationAngle(Real t_tai);
Real GregorianToMJD(int year, int month, int day, int hour = 0, int min = 0,
                    Real sec = 0);
Real GregorianToTime(int year, int month, int day, int hour = 0, int min = 0,
                     Real sec = 0);

Real GreenwichMeanSiderealTime(Real mjd_ut1);
Real GreenwichApparentSiderealTime(Real mjd_ut1);

std::tuple<int, int, int, int, int, Real> MJDtoGregorian(Real mjd);

std::string MJDtoGregorianString(Real mjd, int precision);
std::string TimeToGregorianString(Real t, int precision);

VEC_DEF_REAL(UTCtoUT1)
VEC_DEF_REAL(UT1toUTC)
VEC_DEF_REAL(TAItoUTC)
VEC_DEF_REAL(UTCtoTAI)
VEC_DEF_REAL(TAItoTT)
VEC_DEF_REAL(TTtoTAI)
VEC_DEF_REAL(TCGtoTT)
VEC_DEF_REAL(TTtoTCG)
VEC_DEF_REAL(TTtoTDB)
VEC_DEF_REAL(TDBtoTT)
VEC_DEF_REAL(TAItoGPS)
VEC_DEF_REAL(GPStoTAI)
VEC_DEF_REAL(TCBtoTDB)
VEC_DEF_REAL(TTtoTCB)

VEC_DEF_REAL(MJDtoTime)
VEC_DEF_REAL(TimeToMJD)

}  // namespace lupnt