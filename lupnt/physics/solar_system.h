#pragma once

#include "lupnt/core/constants.h"

namespace lupnt {

Real MeanObliquity(Real mjd_tt);
Mat3 Equatorial2EclipticMatrix(Real mjd_tt);
Mat3 PrecessionMatrix(Real mjd_1, Real mjd_2);
std::pair<Real, Real> NutAngles(Real mjd_tt);
Mat3 NutationMatrix(Real mjd_tt);
Mat3 NutationMatrixLowPrecision(Real mjd_tt);
Real EquinoxEquation(Real mjd_tt);

Vec3 SunPositionLowPrecision(Real mjd_tt);
Vec3 MoonPositionLowPrecision(Real mjd_tt);

Mat3 GreenwichHourAngleMatrix(Real mjd_ut1);

}  // namespace lupnt