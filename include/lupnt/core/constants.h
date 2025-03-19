/**
 * @file constants.h
 * @author Stanford NAV LAB
 * @brief List of constants
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "lupnt/core/definitions.h"
#include "lupnt/core/user_file_path.h"

namespace lupnt {

  // Math constants
  static constexpr double PI = 3.14159265358979323846264338327950288419716939937511;
  static constexpr double TWO_PI = 6.28318530717958647692528676655900576839433879875022;
  static constexpr double PI_OVER_TWO = 1.57079632679489661923132169163975144209858469968756;
  static constexpr double E = 2.71828182845904523536028747135266249775724709369996;
  static constexpr double EPS = 1.0e-16;

  // Angle conversion
  static constexpr double RAD = 3.14159265358979323846264338327950288419716939937511 / 180.0;
  static constexpr double DEG = 180.0 / 3.14159265358979323846264338327950288419716939937511;
  static constexpr double ARCSEC_DEG = 3600.0;
  static constexpr double DEG_ARCSEC = 1.0 / ARCSEC_DEG;
  static constexpr double RAD_ARCSEC = DEG_ARCSEC * RAD;
  static constexpr double ARCSEC_RAD = 1.0 / RAD_ARCSEC;

  // Mass
  static constexpr double LBM_TO_KG = 0.45359237;
  static constexpr double SLUG_TO_KG = 14.59390294;

  // Length
  static constexpr double INCH_M = 0.0254;
  static constexpr double FOOT_M = 0.3048;
  static constexpr double MILE_M = 1609.344;
  static constexpr double KM_M = 0.001;
  static constexpr double M_KM = 1000.0;

  // Time system constants ******************************************************
  static constexpr double SECS_DAY = 86400.0;
  static constexpr double SECS_HOUR = 3600.0;
  static constexpr double SECS_MINUTE = 60.0;

  static constexpr double MINS_HOUR = 60.0;
  static constexpr double MINS_DAY = 1440.0;
  static constexpr double HOURS_DAY = 24.0;
  static constexpr double DAYS_WEEK = 7.0;

  static constexpr double DAYS_YEAR = 365.25;
  static constexpr double DAYS_CENTURY = 36525.00;
  static constexpr double DAYS_SEC = 1.1574074074074074074074074074074e-5;

  static constexpr double TIME_OF_J2000 = 883655990.850000;  // 2000/01/01 43167.85
  static constexpr double JD_J2000 = 2451545.0;              // JD of J2000 epoch
  static constexpr double MJD_J2000 = 51544.5;               // MJD of J2000 epoch

  // Vallado page 94
  static constexpr double JD_T0 = 2443144.5003725;
  static constexpr double JD_MJD_OFFSET = 2400000.5;
  static constexpr double TT_TAI_OFFSET = 32.184;
  static constexpr double A1_TAI_OFFSET = 0.0343817;
  static constexpr double JD_JAN_5_1941 = 2430000.0;
  static constexpr double JD_NOV_17_1858 = 2400000.5;

  static constexpr double L_B = 1.550505e-8;
  static constexpr double L_G = 6.969290134e-10;
  static constexpr double NUM_SECS = SECS_DAY;

  static constexpr int JULIAN_DATE_OF_010541 = 2430000;

  // Coordinate system constants DE440 *******************************************
  static constexpr double GM_SUN = 132712440041.279419;          // [km^3/s^2]
  static constexpr double GM_MERCURY = 22031.868551;             // [km^3/s^2]
  static constexpr double GM_VENUS = 324858.592000;              // [km^3/s^2]
  static constexpr double GM_EARTH = 398600.435507;              // [km^3/s^2]
  static constexpr double GM_MOON = 4902.800118;                 // [km^3/s^2]
  static constexpr double GM_MARS_SYSTEM = 42828.375816;         // [km^3/s^2]
  static constexpr double GM_JUPITER_SYSTEM = 126712764.100000;  // [km^3/s^2]
  static constexpr double GM_SATURN_SYSTEM = 37940584.841800;    // [km^3/s^2]
  static constexpr double GM_URANUS_SYSTEM = 5794556.400000;     // [km^3/s^2]
  static constexpr double GM_NEPTUNE_SYSTEM = 6836527.100580;    // [km^3/s^2]
  static constexpr double GM_PLUTO_SYSTEM = 977.000000;          // [km^3/s^2]
  static constexpr double GM_MARS = 0.4282837566395650E+05;      // [km^3/s^2]
  static constexpr double GM_JUPITER = 0.1267127646799999E+08;   // [km^3/s^2]
  static constexpr double GM_SATURN = 0.3794058480000000E+07;    // [km^3/s^2]
  static constexpr double GM_URANUS = 0.5794556400000000E+06;    // [km^3/s^2]
  static constexpr double GM_NEPTUNE = 0.6836527100580000E+06;   // [km^3/s^2]

  static constexpr double GM_CERES = 62.62890;   // [km^3/s^2]
  static constexpr double GM_VESTA = 17.288245;  // [km^3/s^2]

  static constexpr double D_EARTH_MOON = 384400.0;  // [km]
  static constexpr double D_EARTH_EMB = 4671.0;     // [km]
  static constexpr double R_EARTH = 6378.137;       // [km]
  static constexpr double R_MOON = 1737.4;          // [km]
  static constexpr double R_SUN = 696342.0;         // [km]
  static constexpr double R_MERCURY = 2439.7;       // [km]
  static constexpr double R_VENUS = 6051.8;         // [km]
  static constexpr double R_MARS = 3396.0;          // [km]
  static constexpr double R_JUPITER = 71492.0;      // [km]
  static constexpr double R_SATURN = 60268.0;       // [km]
  static constexpr double R_URANUS = 25559.0;       // [km]
  static constexpr double R_NEPTUNE = 24764.0;      // [km]
  static constexpr double R_PLUTO = 1188.3;         // [km]

  static constexpr double OMEGA_EARTH_MOON = 2.6617e-6;             // [rad/s]
  static constexpr double D_MOON_EMB = D_EARTH_MOON - D_EARTH_EMB;  // [km]

  static constexpr double WGS84_A = 6378.137;             // [km]
  static constexpr double WGS84_F = 1.0 / 298.257223563;  // [-]

  static constexpr double J2_EARTH = 1.08262668e-3;
  // static constexpr double J2_MOON = 9.08901807506000e-5;
  static constexpr double J2_MOON
      = 9.094278450270e-5;  // Zonal value adjusted for permanent tide - Rigid J2
  static constexpr double C22_MOON
      = 3.470983013194e-5;  // Sectorial value adjusted for perm. tide - Rigid C22

  // Transformations Between GCRF and Mean Equator and Equinox at J2000
  static constexpr double FRAME_BIAS_XI0 = -8.0561e-8;     // [rad]
  static constexpr double FRAME_BIAS_ETA0 = -3.3060e-8;    // [rad]
  static constexpr double FRAME_BIAS_DALPHA0 = 7.0783e-8;  // [rad]

  // Solar Radiation Pressure Constants
  static constexpr double AU = 149597970;               // AU [km]
  static constexpr double SOLAR_FLUX_AU = 1358 * 1e-6;  // Mean Solar Flux at 1 AU [W/km^2]
  static constexpr double C = 299792.458;               // Light speed [km/s]
  static constexpr double P_SUN
      = SOLAR_FLUX_AU / C;  // Solar radiation pressure at 1 AU [N/km^2] = 4.56e-6 N/m^2

  // File Path *******************************************************************
  static constexpr std::string_view TAI_UTC_FILENAME = "tai-utc.dat";
  static constexpr std::string_view EOP_FILENAME = "eopc04_08.62-now";
  static constexpr std::string_view IAU_SOFA_FILENAME = "IAU_SOFA.DAT";

  // NAIF Intefer ID codes
  // Reference:
  // https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/naif_ids.html
  enum class NaifId {
    SSB = 0,
    SOLAR_SYSTEM_BARYCENTER = SSB,
    MERCURY_BARYCENTER = 1,
    VENUS_BARYCENTER = 2,
    EMB = 3,
    EARTH_MOON_BARYCENTER = EMB,
    MARS_BARYCENTER = 4,
    JUPITER_BARYCENTER = 5,
    SATURN_BARYCENTER = 6,
    URANUS_BARYCENTER = 7,
    NEPTUNE_BARYCENTER = 8,
    PLUTO_BARYCENTER = 9,
    SUN = 10,
    MERCURY = 199,
    VENUS = 299,
    EARTH = 399,
    MOON = 301,
    MARS = 499,
    PHOBOS = 401,
    DEIMOS = 402,
    JUPITER = 599,
    SATURN = 699,
    URANUS = 799,
    NEPTUNE = 899,
  };

  double GetBodyRadius(NaifId body);
  const std::ostream& operator<<(std::ostream& os, NaifId id);

  enum class Time {
    UT1,    // Universal Time 1
    UTC,    // Coordinated Universal Time
    TAI,    // International Atomic Time
    TDB,    // Barycentric Dynamical Time
    TT,     // Terrestrial Time
    TCG,    // Geocentric Coordinate Time
    TCB,    // Barycentric Coordinate Time
    GPS,    // Global Positioning System Time
    JD_TT,  // Julian Date relative to TT
    JD_TDB  // Julian Date relative to TDB
  };

}  // namespace lupnt
