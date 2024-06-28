#include <string>
#include <tuple>

#include "lupnt/core/constants.h"

namespace lupnt {

Real Calendar2ModJulianDate(int year, int month, int day, int hour = 0,
                            int min = 0, Real sec = 0);

Real GreenwichMeanSiderealTime(Real mjd_ut1);
Real GreenwichApparentSiderealTime(Real mjd_ut1);

std::tuple<int, int, int, int, int, Real> ModJulianDate2Calendar(Real mjd);

std::string FormatDate(Real mjd, int precision);

}  // namespace lupnt