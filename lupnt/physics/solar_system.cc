#include "solar_system.h"

#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/time.h"

namespace lupnt {

/// @brief Compute the mean obliquity of the ecliptic
/// @param mjd_tt Terrestrial Time (Modified Julian Date)
/// @return Mean obliquity of the ecliptic
/// @ref
/// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
/// applications. Berlin : New York: Springer, 2000.
/// doi: 10.1007/978-3-642-58351-3.
Real MeanObliquity(Real mjd_tt) {
  const Real T = (mjd_tt - MJD_J2000) / 36525.0;
  return (23.43929111 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0) *
         RAD;
}

/// @brief Transformation of equatorial to ecliptical coordinates
/// @param mjd_tt Modified Julian Date (Terrestrial Time)
/// @return Transformation matrix
/// @ref
/// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
/// applications. Berlin : New York: Springer, 2000.
/// doi: 10.1007/978-3-642-58351-3.
Mat3 Equatorial2EclipticMatrix(Real mjd_tt) {
  return RotX(MeanObliquity(mjd_tt));
}

/// @brief Precession transformation of equatorial coordinates
/// @param mjd_1 Epoch given (Modified Julian Date TT)
/// @param mjd_2 Epoch to precess to (Modified Julian Date TT)
/// @return Precession transformation matrix
/// @ref
/// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
/// applications. Berlin : New York: Springer, 2000.
/// doi: 10.1007/978-3-642-58351-3.
Mat3 PrecessionMatrix(Real mjd_1, Real mjd_2) {
  Real T = (mjd_1 - MJD_J2000) / 36525.0;
  Real dT = (mjd_2 - mjd_1) / 36525.0;

  // Precession angles
  Real zeta = ((2306.2181 + (1.39656 - 0.000139 * T) * T) +
               ((0.30188 - 0.000344 * T) + 0.017998 * dT) * dT) *
              dT * RAD_ARCSEC;
  Real z =
      zeta + ((0.79280 + 0.000411 * T) + 0.000205 * dT) * dT * dT * RAD_ARCSEC;
  Real theta = ((2004.3109 - (0.85330 + 0.000217 * T) * T) -
                ((0.42665 + 0.000217 * T) + 0.041833 * dT) * dT) *
               dT * RAD_ARCSEC;

  // Precession matrix
  Mat3 P = RotZ(-z) * RotX(theta) * RotZ(-zeta);
  return P;
}

/// @brief Nutation in longitude and obliquity
/// @param mjd_tt Modified Julian Date (Terrestrial Time)
/// @param dpsi Nutation in longitude [rad]
/// @param deps Nutation in obliquity [rad]
/// @return Nutation matrix
/// @ref
/// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
/// applications. Berlin : New York: Springer, 2000.
/// doi: 10.1007/978-3-642-58351-3.
std::pair<Real, Real> NutAngles(Real mjd_tt) {
  const Real T = (mjd_tt - MJD_J2000) / 36525.0;
  const Real T2 = T * T;
  const Real T3 = T2 * T;
  const Real rev = 360.0 * 3600.0;  // arcsec/revolution

  const int N_coeff = 106;
  const long C[N_coeff][9] = {
      //
      // l  l' F  D Om    dpsi    *T     deps     *T       #
      //
      {0, 0, 0, 0, 1, -1719960, -1742, 920250, 89},  //   1
      {0, 0, 0, 0, 2, 20620, 2, -8950, 5},           //   2
      {-2, 0, 2, 0, 1, 460, 0, -240, 0},             //   3
      {2, 0, -2, 0, 0, 110, 0, 0, 0},                //   4
      {-2, 0, 2, 0, 2, -30, 0, 10, 0},               //   5
      {1, -1, 0, -1, 0, -30, 0, 0, 0},               //   6
      {0, -2, 2, -2, 1, -20, 0, 10, 0},              //   7
      {2, 0, -2, 0, 1, 10, 0, 0, 0},                 //   8
      {0, 0, 2, -2, 2, -131870, -16, 57360, -31},    //   9
      {0, 1, 0, 0, 0, 14260, -34, 540, -1},          //  10
      {0, 1, 2, -2, 2, -5170, 12, 2240, -6},         //  11
      {0, -1, 2, -2, 2, 2170, -5, -950, 3},          //  12
      {0, 0, 2, -2, 1, 1290, 1, -700, 0},            //  13
      {2, 0, 0, -2, 0, 480, 0, 10, 0},               //  14
      {0, 0, 2, -2, 0, -220, 0, 0, 0},               //  15
      {0, 2, 0, 0, 0, 170, -1, 0, 0},                //  16
      {0, 1, 0, 0, 1, -150, 0, 90, 0},               //  17
      {0, 2, 2, -2, 2, -160, 1, 70, 0},              //  18
      {0, -1, 0, 0, 1, -120, 0, 60, 0},              //  19
      {-2, 0, 0, 2, 1, -60, 0, 30, 0},               //  20
      {0, -1, 2, -2, 1, -50, 0, 30, 0},              //  21
      {2, 0, 0, -2, 1, 40, 0, -20, 0},               //  22
      {0, 1, 2, -2, 1, 40, 0, -20, 0},               //  23
      {1, 0, 0, -1, 0, -40, 0, 0, 0},                //  24
      {2, 1, 0, -2, 0, 10, 0, 0, 0},                 //  25
      {0, 0, -2, 2, 1, 10, 0, 0, 0},                 //  26
      {0, 1, -2, 2, 0, -10, 0, 0, 0},                //  27
      {0, 1, 0, 0, 2, 10, 0, 0, 0},                  //  28
      {-1, 0, 0, 1, 1, 10, 0, 0, 0},                 //  29
      {0, 1, 2, -2, 0, -10, 0, 0, 0},                //  30
      {0, 0, 2, 0, 2, -22740, -2, 9770, -5},         //  31
      {1, 0, 0, 0, 0, 7120, 1, -70, 0},              //  32
      {0, 0, 2, 0, 1, -3860, -4, 2000, 0},           //  33
      {1, 0, 2, 0, 2, -3010, 0, 1290, -1},           //  34
      {1, 0, 0, -2, 0, -1580, 0, -10, 0},            //  35
      {-1, 0, 2, 0, 2, 1230, 0, -530, 0},            //  36
      {0, 0, 0, 2, 0, 630, 0, -20, 0},               //  37
      {1, 0, 0, 0, 1, 630, 1, -330, 0},              //  38
      {-1, 0, 0, 0, 1, -580, -1, 320, 0},            //  39
      {-1, 0, 2, 2, 2, -590, 0, 260, 0},             //  40
      {1, 0, 2, 0, 1, -510, 0, 270, 0},              //  41
      {0, 0, 2, 2, 2, -380, 0, 160, 0},              //  42
      {2, 0, 0, 0, 0, 290, 0, -10, 0},               //  43
      {1, 0, 2, -2, 2, 290, 0, -120, 0},             //  44
      {2, 0, 2, 0, 2, -310, 0, 130, 0},              //  45
      {0, 0, 2, 0, 0, 260, 0, -10, 0},               //  46
      {-1, 0, 2, 0, 1, 210, 0, -100, 0},             //  47
      {-1, 0, 0, 2, 1, 160, 0, -80, 0},              //  48
      {1, 0, 0, -2, 1, -130, 0, 70, 0},              //  49
      {-1, 0, 2, 2, 1, -100, 0, 50, 0},              //  50
      {1, 1, 0, -2, 0, -70, 0, 0, 0},                //  51
      {0, 1, 2, 0, 2, 70, 0, -30, 0},                //  52
      {0, -1, 2, 0, 2, -70, 0, 30, 0},               //  53
      {1, 0, 2, 2, 2, -80, 0, 30, 0},                //  54
      {1, 0, 0, 2, 0, 60, 0, 0, 0},                  //  55
      {2, 0, 2, -2, 2, 60, 0, -30, 0},               //  56
      {0, 0, 0, 2, 1, -60, 0, 30, 0},                //  57
      {0, 0, 2, 2, 1, -70, 0, 30, 0},                //  58
      {1, 0, 2, -2, 1, 60, 0, -30, 0},               //  59
      {0, 0, 0, -2, 1, -50, 0, 30, 0},               //  60
      {1, -1, 0, 0, 0, 50, 0, 0, 0},                 //  61
      {2, 0, 2, 0, 1, -50, 0, 30, 0},                //  62
      {0, 1, 0, -2, 0, -40, 0, 0, 0},                //  63
      {1, 0, -2, 0, 0, 40, 0, 0, 0},                 //  64
      {0, 0, 0, 1, 0, -40, 0, 0, 0},                 //  65
      {1, 1, 0, 0, 0, -30, 0, 0, 0},                 //  66
      {1, 0, 2, 0, 0, 30, 0, 0, 0},                  //  67
      {1, -1, 2, 0, 2, -30, 0, 10, 0},               //  68
      {-1, -1, 2, 2, 2, -30, 0, 10, 0},              //  69
      {-2, 0, 0, 0, 1, -20, 0, 10, 0},               //  70
      {3, 0, 2, 0, 2, -30, 0, 10, 0},                //  71
      {0, -1, 2, 2, 2, -30, 0, 10, 0},               //  72
      {1, 1, 2, 0, 2, 20, 0, -10, 0},                //  73
      {-1, 0, 2, -2, 1, -20, 0, 10, 0},              //  74
      {2, 0, 0, 0, 1, 20, 0, -10, 0},                //  75
      {1, 0, 0, 0, 2, -20, 0, 10, 0},                //  76
      {3, 0, 0, 0, 0, 20, 0, 0, 0},                  //  77
      {0, 0, 2, 1, 2, 20, 0, -10, 0},                //  78
      {-1, 0, 0, 0, 2, 10, 0, -10, 0},               //  79
      {1, 0, 0, -4, 0, -10, 0, 0, 0},                //  80
      {-2, 0, 2, 2, 2, 10, 0, -10, 0},               //  81
      {-1, 0, 2, 4, 2, -20, 0, 10, 0},               //  82
      {2, 0, 0, -4, 0, -10, 0, 0, 0},                //  83
      {1, 1, 2, -2, 2, 10, 0, -10, 0},               //  84
      {1, 0, 2, 2, 1, -10, 0, 10, 0},                //  85
      {-2, 0, 2, 4, 2, -10, 0, 10, 0},               //  86
      {-1, 0, 4, 0, 2, 10, 0, 0, 0},                 //  87
      {1, -1, 0, -2, 0, 10, 0, 0, 0},                //  88
      {2, 0, 2, -2, 1, 10, 0, -10, 0},               //  89
      {2, 0, 2, 2, 2, -10, 0, 0, 0},                 //  90
      {1, 0, 0, 2, 1, -10, 0, 0, 0},                 //  91
      {0, 0, 4, -2, 2, 10, 0, 0, 0},                 //  92
      {3, 0, 2, -2, 2, 10, 0, 0, 0},                 //  93
      {1, 0, 2, -2, 0, -10, 0, 0, 0},                //  94
      {0, 1, 2, 0, 1, 10, 0, 0, 0},                  //  95
      {-1, -1, 0, 2, 1, 10, 0, 0, 0},                //  96
      {0, 0, -2, 0, 1, -10, 0, 0, 0},                //  97
      {0, 0, 2, -1, 2, -10, 0, 0, 0},                //  98
      {0, 1, 0, 2, 0, -10, 0, 0, 0},                 //  99
      {1, 0, -2, -2, 0, -10, 0, 0, 0},               // 100
      {0, -1, 2, 0, 1, -10, 0, 0, 0},                // 101
      {1, 1, 0, -2, 1, -10, 0, 0, 0},                // 102
      {1, 0, -2, 2, 0, -10, 0, 0, 0},                // 103
      {2, 0, 0, 2, 0, 10, 0, 0, 0},                  // 104
      {0, 0, 2, 4, 2, -10, 0, 0, 0},                 // 105
      {0, 1, 0, 1, 0, 10, 0, 0, 0}                   // 106
  };

  // Mean arguments of luni-solar motion
  //
  //   l   mean anomaly of the Moon
  //   l'  mean anomaly of the Sun
  //   F   mean argument of latitude
  //   D   mean longitude elongation of the Moon from the Sun
  //   Om  mean longitude of the ascending node
  Real l = mod(
      485866.733 + (1325.0 * rev + 715922.633) * T + 31.310 * T2 + 0.064 * T3,
      rev);
  Real lp = mod(
      1287099.804 + (99.0 * rev + 1292581.224) * T - 0.577 * T2 - 0.012 * T3,
      rev);
  Real F = mod(
      335778.877 + (1342.0 * rev + 295263.137) * T - 13.257 * T2 + 0.011 * T3,
      rev);
  Real D = mod(
      1072261.307 + (1236.0 * rev + 1105601.328) * T - 6.891 * T2 + 0.019 * T3,
      rev);
  Real Om = mod(
      450160.280 - (5.0 * rev + 482890.539) * T + 7.455 * T2 + 0.008 * T3, rev);

  // Nutation in longitude and obliquity [rad]
  Real arg;
  Real deps = 0, dpsi = 0;
  for (int i = 0; i < N_coeff; i++) {
    arg = (C[i][0] * l + C[i][1] * lp + C[i][2] * F + C[i][3] * D +
           C[i][4] * Om) *
          RAD_ARCSEC;
    dpsi += (C[i][5] + C[i][6] * T) * sin(arg);
    deps += (C[i][7] + C[i][8] * T) * cos(arg);
  };

  dpsi = 1.0E-5 * dpsi * RAD_ARCSEC;
  deps = 1.0E-5 * deps * RAD_ARCSEC;
  return {dpsi, deps};
}

/// @brief Transformation from mean to true equator and equinox
/// @param mjd_tt Terrestrial Time (Modified Julian Date)
/// @return Nutation matrix
/// @ref
/// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
/// applications. Berlin : New York: Springer, 2000.
/// doi: 10.1007/978-3-642-58351-3.
Mat3 NutationMatrix(Real mjd_tt) {
  // Mean obliquity of the ecliptic
  Real eps = MeanObliquity(mjd_tt);

  // Nutation in longitude and obliquity
  auto [dpsi, deps] = NutAngles(mjd_tt);

  // Transformation from mean to true equator and equinox
  Mat3 R = RotX(-eps - deps) * RotZ(-dpsi) * RotX(+eps);
  return R;
}

/// @brief Transformation from mean to true equator and equinox (low precision)
/// @param mjd_tt Terrestrial Time (Modified Julian Date)
/// @return Nutation matrix
/// @ref
/// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
/// applications. Berlin : New York: Springer, 2000.
/// doi: 10.1007/978-3-642-58351-3.
Mat3 NutationMatrixLowPrecision(Real mjd_tt) {
  Real T = (mjd_tt - MJD_J2000) / 36525.0;

  // Mean arguments of luni-solar motion
  //   ls  mean anomaly of the Sun
  //   D   diff. longitude Moon-Sun
  //   F   mean argument of latitude
  //   N   longit. ascending node
  Real ls = TWO_PI * frac(0.993133 + 99.997306 * T);
  Real D = TWO_PI * frac(0.827362 + 1236.853087 * T);
  Real F = TWO_PI * frac(0.259089 + 1342.227826 * T);
  Real N = TWO_PI * frac(0.347346 - 5.372447 * T);

  // Nutation angles
  Real dpsi =
      (-17.200 * sin(N) - 1.319 * sin(2 * (F - D + N)) -
       0.227 * sin(2 * (F + N)) + 0.206 * sin(2 * N) + 0.143 * sin(ls)) *
      RAD_ARCSEC;
  Real deps = (+9.203 * cos(N) + 0.574 * cos(2 * (F - D + N)) +
               0.098 * cos(2 * (F + N)) - 0.090 * cos(2 * N)) *
              RAD_ARCSEC;

  // Mean obliquity of the ecliptic
  Real eps = 0.4090928 - 2.2696E-4 * T;
  Mat3 R = RotX(-eps - deps) * RotZ(-dpsi) * RotX(+eps);
  return R;
}

/// @brief Computation of the equation of the equinoxes
/// @param mjd_tt Terrestrial Time (Modified Julian Date)
/// @return Equation of the equinoxes
/// @note The equation of the equinoxes dpsi*cos(eps) is the right ascension of
/// the mean equinox referred to the true equator and equinox and is equal to
/// the difference between apparent and mean sidereal time.
/// @ref
/// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
/// applications. Berlin : New York: Springer, 2000.
/// doi: 10.1007/978-3-642-58351-3.
Real EquinoxEquation(Real mjd_tt) {
  // Nutation in longitude and obliquity
  auto [dpsi, reps] = NutAngles(mjd_tt);
  // Equation of the equinoxes
  return dpsi * cos(MeanObliquity(mjd_tt));
};

/// @brief Compute the Sun's geocentric position using a low precision
/// @param mjd_tt Terrestrial Time (Modified Julian Date)
/// @return Solar position vector [km] with respect to the mean equator and
/// equinox of J2000 (EME2000, ICRF)
/// @ref
/// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
/// applications. Berlin : New York: Springer, 2000.
/// doi: 10.1007/978-3-642-58351-3.
Vec3 SunPositionLowPrecision(Real mjd_tt) {
  const Real eps = 23.43929111 * RAD;             // Obliquity of J2000 ecliptic
  const Real T = (mjd_tt - MJD_J2000) / 36525.0;  // Julian cent. since J2000

  // Mean anomaly, ecliptic longitude and radius
  Real M = TWO_PI * frac(0.9931267 + 99.9973583 * T);  // [rad]
  Real L = TWO_PI *
           frac(0.7859444 + M / TWO_PI +
                (6892.0 * sin(M) + 72.0 * sin(2.0 * M)) / 1296.0e3);  // [rad]
  Real r = 149.619e9 - 2.499e9 * cos(M) - 0.021e9 * cos(2 * M);       // [m]

  // Equatorial position vector
  Vec3 r_sun = RotX(-eps) * Vec3(r * cos(L), r * sin(L), 0.0);
  return r_sun;
}

/// @brief Compute the Moon's geocentric position using a low precision
/// analytical series
/// @param mjd_tt Terrestrial Time (Modified Julian Date)
/// @return Lunar position vector [km] with respect to the mean equator and
/// equinox of J2000 (EME2000, ICRF)
/// @ref
/// O. Montenbruck and G. Eberhard, Satellite orbits: models, methods, and
/// applications. Berlin : New York: Springer, 2000.
/// doi: 10.1007/978-3-642-58351-3.
Vec3 MoonPositionLowPrecision(Real mjd_tt) {
  const Real eps = 23.43929111 * RAD;             // Obliquity of J2000 ecliptic
  const Real T = (mjd_tt - MJD_J2000) / 36525.0;  // Julian cent. since J2000

  // Mean elements of lunar orbit
  Real L_0 = frac(
      0.606433 + 1336.851344 * T);  // Mean longitude [rev] w.r.t. J2000 equinox
  Real l =
      TWO_PI * frac(0.374897 + 1325.552410 * T);  // Moon's mean anomaly [rad]
  Real lp =
      TWO_PI * frac(0.993133 + 99.997361 * T);  // Sun's mean anomaly [rad]
  Real D =
      TWO_PI * frac(0.827361 + 1236.853086 * T);  // Diff. long. Moon-Sun [rad]
  Real F = TWO_PI * frac(0.259086 + 1342.227825 * T);  // Argument of latitude

  // Ecliptic longitude (w.r.t. equinox of J2000)
  Real dL = +22640 * sin(l) - 4586 * sin(l - 2 * D) + 2370 * sin(2 * D) +
            769 * sin(2 * l) - 668 * sin(lp) - 412 * sin(2 * F) -
            212 * sin(2 * l - 2 * D) - 206 * sin(l + lp - 2 * D) +
            192 * sin(l + 2 * D) - 165 * sin(lp - 2 * D) - 125 * sin(D) -
            110 * sin(l + lp) + 148 * sin(l - lp) - 55 * sin(2 * F - 2 * D);
  Real L = TWO_PI * frac(L_0 + dL / 1296.0e3);  // [rad]

  // Ecliptic latitude
  Real S = F + (dL + 412 * sin(2 * F) + 541 * sin(lp)) * RAD_ARCSEC;
  Real h = F - 2 * D;
  Real N = -526 * sin(h) + 44 * sin(l + h) - 31 * sin(-l + h) -
           23 * sin(lp + h) + 11 * sin(-lp + h) - 25 * sin(-2 * l + F) +
           21 * sin(-l + F);
  Real B = (18520.0 * sin(S) + N) * RAD_ARCSEC;  // [rad]
  Real cosB = cos(B);

  // Distance [m]
  Real R = 385000e3 - 20905e3 * cos(l) - 3699e3 * cos(2 * D - l) -
           2956e3 * cos(2 * D) - 570e3 * cos(2 * l) +
           246e3 * cos(2 * l - 2 * D) - 205e3 * cos(lp - 2 * D) -
           171e3 * cos(l + 2 * D) - 152e3 * cos(l + lp - 2 * D);

  // Equatorial coordinates
  Vec3 r_moon =
      RotX(-eps) * Vec3(R * cos(L) * cosB, R * sin(L) * cosB, R * sin(B));
  return r_moon;
}

//------------------------------------------------------------------------------
//
// GHAMatrix
//
// Purpose:
//
//   Transformation from true equator and equinox to Earth equator and
//   Greenwich meridian system
//
// Input/Output:
//
//   mjd_ut1   Modified Julian Date UT1
//   <return>  Greenwich Hour Angle matrix
//
//------------------------------------------------------------------------------

/// @brief Transformation from true equator and equinox to Earth equator and
/// Greenwich meridian system
/// @param mjd_ut1 UT1 (Modified Julian Date)
/// @return Greenwich Hour Angle matrix
Mat3 GreenwichHourAngleMatrix(Real mjd_ut1) {
  return RotZ(GreenwichApparentSiderealTime(mjd_ut1));
}

}  // namespace lupnt