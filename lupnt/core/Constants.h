#pragma once

#include <autodiff/common/eigen.hpp>
#include <autodiff/forward/real.hpp>
#include <string>

#include "UserFilePath.h"

namespace Eigen {
typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;

}  // namespace Eigen

namespace autodiff {
using Matrix6real = Eigen::Matrix<real1st, 6, 6, 0, 6, 6>;
using Vector6real = Eigen::Matrix<real1st, 6, 1, 0, 6, 1>;
}  // namespace autodiff

namespace LPT {

// Math constants --------------------------------------------------------------
static constexpr double PI_DEG = 180.0;
static constexpr double PI_OVER_TWO_DEG = 90.0;
static constexpr double TWO_PI_DEG = 360.0;
static constexpr double PI =
    3.14159265358979323846264338327950288419716939937511;
static constexpr double TWO_PI =
    6.28318530717958647692528676655900576839433879875022;
static constexpr double PI_OVER_TWO =
    1.57079632679489661923132169163975144209858469968756;
static constexpr double E =
    2.71828182845904523536028747135266249775724709369996;

// Angle conversion
static constexpr double RAD_PER_DEG =
    3.14159265358979323846264338327950288419716939937511 / 180.0;
static constexpr double DEG_PER_RAD =
    180.0 / 3.14159265358979323846264338327950288419716939937511;

static constexpr double ARCSEC_PER_DEGREE = 3600.0;
static constexpr double DEG_PER_ARCSEC = 1.0 / 3600.0;
static constexpr double RAD_PER_ARCSEC = DEG_PER_ARCSEC * RAD_PER_DEG;

// Mass conversion
static constexpr double LBM_TO_KG = 0.45359237;
static constexpr double SLUG_TO_KG = 14.59390294;

// Length
static constexpr double INCH_TO_M = 0.0254;
static constexpr double FOOT_TO_M = 0.3048;
static constexpr double STATUTE_MILE_TO_M = 1609.344;
static constexpr double M_TO_KM = 0.001;
static constexpr double KM_TO_M = 1000.0;

// Time system constants -------------------------------------------------------
static constexpr double SECS_PER_DAY = 86400.0;
static constexpr double SECS_PER_HOUR = 3600.0;
static constexpr double SECS_PER_MINUTE = 60.0;

static constexpr double DAYS_PER_YEAR = 365.25;
static constexpr double DAYS_PER_JULIAN_CENTURY = 36525.00;
static constexpr double DAYS_PER_SEC = 1.1574074074074074074074074074074e-5;

static constexpr double TIME_OF_J2000 =
    883655990.850000;                                   // 2000/01/01 43167.85
static constexpr double JD_OF_J2000 = 2451545.0;        // JD of J2000 epoch
static constexpr double MJD_OF_J2000 = 21545.00000000;  // MJD of J2000 epoch
static constexpr double A1MJD_OF_J2000 =
    21545.00000000;  // 2000/01/01 11:59:27.965622
static constexpr double JD_MJD_OFFSET =
    2400000.5;  // Vallado page 187 (= JD_NOV_17_1858)
static constexpr double TT_TAI_OFFSET = 32.184;  // GMAT Math Spec section 2.3
static constexpr double A1_TAI_OFFSET =
    0.0343817;  // GMAT Math Spec section 2.1
static constexpr double JD_JAN_5_1941 =
    2430000.0;  // old name JULIAN_DATE_OF_010541
static constexpr double JD_NOV_17_1858 = 2400000.5;  // old name JD_MJD_OFFSET

static constexpr double TDB_COEFF1 = 0.001658;
static constexpr double TDB_COEFF2 = 0.00001385;
static constexpr double M_E_OFFSET = 357.5277233;
static constexpr double M_E_COEFF1 = 35999.05034;
static constexpr double T_TT_OFFSET = JD_OF_J2000;
static constexpr double T_TT_COEFF1 = DAYS_PER_JULIAN_CENTURY;
static constexpr double L_B = 1.550505e-8;
static constexpr double NUM_SECS = SECS_PER_DAY;

static constexpr int JULIAN_DATE_OF_010541 = 2430000;

// Coordinate system constants -------------------------------------------------
static constexpr double d_E_M = 384400.0;        // km
static constexpr double MU_EARTH = 398600.4418;  // km^3/s^2
static constexpr double MU_MOON = 4902.800066;   // km^3/s^2
static constexpr double d_E_EMB = 4671.0;        // km
static constexpr double R_EARTH = 6378.137;      // km
static constexpr double R_MOON = 1737.4;         // km
static constexpr double OMEGA_E_M = 2.6617e-6;   // rad/s
static constexpr double d_M_EMB = d_E_M - d_E_EMB;

static constexpr double J2_EARTH = 1.08262668e-3;
// static constexpr double J2_MOON = 9.08901807506000e-5;
static constexpr double J2_MOON =
    9.094278450270e-5;  // Zonal value adjusted for permanent tide - Rigid J2
static constexpr double C22_MOON =
    3.470983013194e-5;  // Sectorial value adjusted for perm. tide - Rigid C22

// Solar Radiation Pressure Constants
// -------------------------------------------------
static constexpr double AU = 149597970;      // AU [km]
static constexpr double S_AU = 1358 * 1e-6;  // Mean Solar Flux at 1 AU [W/km^2]
static constexpr double C = 299792.458;      // Light speed [km/s]
static constexpr double P_SUN =
    S_AU / C;  // Solar radiation pressure at 1 AU [N/km^2] = 4.56e-6 N/m^2

// File Pathes -----------------------------------------------------------------
static const std::filesystem::path CSPICE_KER_DIR =
    BASEPATH / "data" / "ephemeris";

// Moon mean elements

}  // namespace LPT