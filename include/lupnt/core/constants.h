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

#include "user_file_path.h"

#define ASSERT_WITH_MESSAGE(condition, message) \
  if (!(condition)) {                           \
    std::ostringstream oss;                     \
    oss << message;                             \
    throw std::runtime_error(oss.str());        \
  }

#define DEFINE_STATIC_VECTORS_MATRICES(size)       \
  using Vec##size = Matrix<Real, size, 1>;         \
  using Vec##size##d = Matrix<double, size, 1>;    \
  using Vec##size##i = Matrix<int, size, 1>;       \
  using Mat##size = Matrix<Real, size, size>;      \
  using Mat##size##d = Matrix<double, size, size>; \
  using Mat##size##i = Matrix<int, size, size>;    \
  using RowVec##size##d = Matrix<double, 1, size>; \
  using RowVec##size##i = Matrix<int, 1, size>;    \
  using RowVec##size = Matrix<Real, 1, size>;

#define DEFINE_DYNAMIC_VECTORS_MATRICES()                       \
  using VecX = Matrix<Real, Eigen::Dynamic, 1>;                 \
  using VecXd = Matrix<double, Eigen::Dynamic, 1>;              \
  using VecXi = Matrix<int, Eigen::Dynamic, 1>;                 \
  using MatX = Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;    \
  using MatXd = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; \
  using MatXi = Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;    \
  using RowVecX = Matrix<Real, 1, Eigen::Dynamic>;              \
  using RowVecXd = Matrix<double, 1, Eigen::Dynamic>;           \
  using RowVecXi = Matrix<int, 1, Eigen::Dynamic>;

#define DEFINE_VECTORS_MATRICES()    \
  DEFINE_STATIC_VECTORS_MATRICES(1)  \
  DEFINE_STATIC_VECTORS_MATRICES(2)  \
  DEFINE_STATIC_VECTORS_MATRICES(3)  \
  DEFINE_STATIC_VECTORS_MATRICES(4)  \
  DEFINE_STATIC_VECTORS_MATRICES(5)  \
  DEFINE_STATIC_VECTORS_MATRICES(6)  \
  DEFINE_STATIC_VECTORS_MATRICES(7)  \
  DEFINE_STATIC_VECTORS_MATRICES(8)  \
  DEFINE_STATIC_VECTORS_MATRICES(9)  \
  DEFINE_STATIC_VECTORS_MATRICES(10) \
  DEFINE_DYNAMIC_VECTORS_MATRICES()

namespace lupnt {

  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixX;
  using Eigen::Vector;
  using Eigen::VectorX;

  using Real = autodiff::real;
  template <int rows, int cols> using Mat = Matrix<Real, rows, cols>;
  template <int rows, int cols> using Matd = Matrix<double, rows, cols>;
  template <int size> using Vec = Matrix<Real, size, 1>;
  template <int size> using Vecd = Matrix<double, size, 1>;
  template <int size> using RowVecd = Matrix<double, 1, size>;
  using RowVecXd = Matrix<double, 1, Eigen::Dynamic>;
  using Quat = Eigen::Quaternion<Real>;
  using Quatd = Eigen::Quaternion<double>;
  using AngleAxis = Eigen::AngleAxis<Real>;
  using AngleAxisd = Eigen::AngleAxis<double>;

  DEFINE_VECTORS_MATRICES()

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
  // Vallado page 187 (= JD_NOV_17_1858)
  static constexpr double JD_MJD_OFFSET = 2400000.5;
  // GMAT Math Spec section 2.3
  static constexpr double TT_TAI_OFFSET = 32.184;
  // GMAT Math Spec section 2.1
  static constexpr double A1_TAI_OFFSET = 0.0343817;
  // old name JULIAN_DATE_OF_010541
  static constexpr double JD_JAN_5_1941 = 2430000.0;
  // old name JD_MJD_OFFSET
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
  static constexpr double GM_PLUTO_SYSTEM = 975.500000;          // [km^3/s^2]
  static constexpr double GM_CERES = 62.62890;                   // [km^3/s^2]
  static constexpr double GM_VESTA = 17.288245;                  // [km^3/s^2]

  static constexpr double D_EARTH_MOON = 384400.0;                  // [km]
  static constexpr double D_EARTH_EMB = 4671.0;                     // [km]
  static constexpr double R_EARTH = 6378.137;                       // [km]
  static constexpr double R_MOON = 1737.4;                          // [km]
  static constexpr double OMEGA_E_M = 2.6617e-6;                    // [rad/s]
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
  static constexpr double AU = 149597970;      // AU [km]
  static constexpr double S_AU = 1358 * 1e-6;  // Mean Solar Flux at 1 AU [W/km^2]
  static constexpr double C = 299792.458;      // Light speed [km/s]
  static constexpr double P_SUN
      = S_AU / C;  // Solar radiation pressure at 1 AU [N/km^2] = 4.56e-6 N/m^2

  // File Path *******************************************************************
  static const std::filesystem::path CSPICE_KER_DIR = GetDataPath() / "ephemeris";
  constexpr const char* TAI_UTC_FILENAME = "tai-utc.dat";
  constexpr const char* EOP_FILENAME = "eopc04_08.62-now";
  constexpr const char* IAU_SOFA_FILENAME = "IAU_SOFA.DAT";
  static const std::filesystem::path ASCII_KERNEL_DIR = GetDataPath() / "ephemeris" / "ascii";

  // Moon mean elements

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

  static std::string toString(NaifId id) {
    switch (id) {
      case NaifId::SOLAR_SYSTEM_BARYCENTER:
        return "SOLAR_SYSTEM_BARYCENTER";
      case NaifId::MERCURY_BARYCENTER:
        return "MERCURY_BARYCENTER";
      case NaifId::VENUS_BARYCENTER:
        return "VENUS_BARYCENTER";
      case NaifId::EARTH_MOON_BARYCENTER:
        return "EARTH_MOON_BARYCENTER";
      case NaifId::MARS_BARYCENTER:
        return "MARS_BARYCENTER";
      case NaifId::JUPITER_BARYCENTER:
        return "JUPITER_BARYCENTER";
      case NaifId::SATURN_BARYCENTER:
        return "SATURN_BARYCENTER";
      case NaifId::URANUS_BARYCENTER:
        return "URANUS_BARYCENTER";
      case NaifId::NEPTUNE_BARYCENTER:
        return "NEPTUNE_BARYCENTER";
      case NaifId::PLUTO_BARYCENTER:
        return "PLUTO_BARYCENTER";
      case NaifId::SUN:
        return "SUN";
      case NaifId::MERCURY:
        return "MERCURY";
      case NaifId::VENUS:
        return "VENUS";
      case NaifId::EARTH:
        return "EARTH";
      case NaifId::MOON:
        return "MOON";
      case NaifId::MARS:
        return "MARS";
      case NaifId::PHOBOS:
        return "PHOBOS";
      case NaifId::DEIMOS:
        return "DEIMOS";
      case NaifId::JUPITER:
        return "JUPITER";
      case NaifId::SATURN:
        return "SATURN";
      case NaifId::URANUS:
        return "URANUS";
      case NaifId::NEPTUNE:
        return "NEPTUNE";
      default:
        return "UNKNOWN";
    }
  }

  namespace TimeSys {
    constexpr const char* UT1 = "UT1";  // Universal Time 1
    constexpr const char* UTC = "UTC";  // Coordinated Universal Time
    constexpr const char* TAI = "TAI";  // International Atomic Time
    constexpr const char* TDB = "TDB";  // Barycentric Dynamical Time
    constexpr const char* TT = "TT";    // Terrestrial Time
    constexpr const char* TCG = "TCG";  // Geocentric Coordinate Time
    constexpr const char* TCB = "TCB";  // Barycentric Coordinate Time
    constexpr const char* GPS = "GPS";  // Global Positioning System Time
    constexpr const char* JD_TT = "JDTDT";
    constexpr const char* JD_TDB = "JDTDB";
  }  // namespace TimeSys

}  // namespace lupnt