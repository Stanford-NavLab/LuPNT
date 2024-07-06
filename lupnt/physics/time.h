#include <string>
#include <tuple>

#include "lupnt/core/constants.h"

namespace lupnt {

Real EarthRotationAngle(Real mjd_ut1);

Real TAItoTT(Real tai);
Real TTtoTAI(Real tt);

Real TTtoTCG(Real tt);
Real TCGtoTT(Real tcg);

Real UTCtoUT1(Real mjd_utc);

Real TAItoJulianDateTT(Real tai);
Real TTtoTDB(Real tt, Real jdtt);
Real TAItoTDB(Real tai);

Real Calendar2ModJulianDate(int year, int month, int day, int hour = 0,
                            int min = 0, Real sec = 0);

Real GreenwichMeanSiderealTime(Real mjd_ut1);
Real GreenwichApparentSiderealTime(Real mjd_ut1);

std::tuple<int, int, int, int, int, Real> ModJulianDate2Calendar(Real mjd);

std::string FormatDate(Real mjd, int precision);

}  // namespace lupnt