/**
 * @file conversions.cpp
 * @author Keidai Iiyama
 * @brief  This file contains the implementation of the coordinate conversion
 * functions
 * @version 0.1
 * @date 2025-02-07
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "lupnt/environment/plasma/gcpm/conversions.h"

#include "lupnt/environment/plasma/core/math_utils.h"
#include "lupnt/environment/plasma/gcpm/constants_gcpm.h"

namespace pecsim {
  void pol_to_cart(double lat, double lon, double r, std::array<double, 3>& pos_sm) {
    double clat = cos(lat);
    pos_sm[0] = r * clat * cos(lon);
    pos_sm[1] = r * clat * sin(lon);
    pos_sm[2] = r * sin(lat);
  }

  void cart_to_pol(const std::array<double, 3>& pos_geo, double& blatr, double& blong,
                   double& rtemp) {
    double x = pos_geo[0];
    double y = pos_geo[1];
    double z = pos_geo[2];
    double row = sqrt(x * x + y * y);
    rtemp = sqrt(x * x + y * y + z * z);
    blatr = atan2(z, row);
    blong = atan2(y, x);
    blong = blong + (1.0 - sign(1.0, blong)) * PI;
  }

  /**
   * Frames:
   * GEI: Geocentric Equatorial Inertial: An inertial (non‐rotating) reference
   * frame centered at the Earth. x: points to the vernal equinox  z: aligned with
   * the Earth’s rotation axis GEO: Geocentric Earth-Fixed Frame: A rotating
   * reference frame centered at the Earth. (Geographic coordinates) x: points
   * towards intersection of the prime meridian and the equator  z: aligned with
   * the Earth’s rotation axis GSE: Geocentric Solar Ecliptic Frame: A rotating
   * reference frame centered at the Earth. x: points towards the Sun  z:
   * perpendicular to the ecliptic plane ("north" of the Earth orbit) GSM:
   * Geocentric Solar Magnetospheric Frame: A rotating reference frame centered at
   * the Earth. x: points towards the Sun  z: aligned with the magnetic dipole
   * axis SM: Solar Magnetic Frame: A rotating reference frame centered at the
   * Earth. x: plane defined by magnetic dipole axis and the Earth-Sun line  z:
   * aligned with the magnetic dipole axis
   */

  // transformation from SM to GEO coordinates
  // SM coordinates are the coordinates defined in the Solar Magnetic (SM)
  // coordinate system and used in the GCPM model Geographic corrdinates are the
  // coordinates defined in ECEF and are inputs for the IRI model
  void sm_to_geo(const std::array<int, 2>& itime, const std::array<double, 3>& x_in,
                 std::array<double, 3>& x_out) {
    std::array<double, 3> temp1, temp2, temp3;

    t4(itime, x_in, temp1, -1);
    t3(itime, temp1, temp2, -1);
    t2(itime, temp2, temp3, -1);
    t1(itime, temp3, x_out, 1);
  }

  void geo_to_sm(const std::array<int, 2>& itime, const std::array<double, 3>& x_in,
                 std::array<double, 3>& x_out) {
    std::array<double, 3> temp1, temp2, temp3;

    t1(itime, x_in, temp1, -1);  // GEO -> GEI
    t2(itime, temp1, temp2, 1);  // GEI -> GSE
    t3(itime, temp2, temp3, 1);  // GSE -> GSM
    t4(itime, temp3, x_out, 1);  // GSM -> SM
  }

  void gei_to_sm(const std::array<int, 2>& itime, const std::array<double, 3>& x_in,
                 std::array<double, 3>& x_out) {
    std::array<double, 3> temp1, temp2;

    t2(itime, x_in, temp1, 1);   // GEI -> GSE
    t3(itime, temp1, temp2, 1);  // GSE -> GSM
    t4(itime, temp2, x_out, 1);  // GSM -> SM
  }

  // function subroutine for computing time in Julian centuries
  double to_mjd(const std::array<int, 2>& itime, double& iyr, double& iday, double& ut) {
    iyr = itime[0] / 1000;
    iday = itime[0] - iyr * 1000;
    ut = itime[1] / 1000.0 / 3600.0;  // UT in hours
    double fracday = ut / 24.0;

    double rmjd = 45.0 + double(iyr - 1859) * 365.0 + (float((iyr - 1861) / 4) + 1.0) + double(iday)
                  - 1.0 + fracday;
    double T0 = (rmjd - 51544.5) / 36525.0;

    return T0;
  }

  // transformation matrix t1=<theta,z>
  // GEI -> GEO
  // Rotates by the Greenwich sidereal time (θ) about the z‑axis to “lock” the
  // inertial frame to the rotating Earth.
  void t1(const std::array<int, 2>& itime, const std::array<double, 3>& x_in,
          std::array<double, 3>& x_out, int iverse) {
    double iyr, iday, ut;
    double temp = to_mjd(itime, iyr, iday, ut);
    double theta = (100.461 + 36000.770 * temp + 15.04107 * ut) * DEG2RAD;
    rotate_z(double(iverse) * theta, x_in, x_out);
  }

  // transformation matrix t2=<lamdas,Z>*<epsilon,X>
  // GEI -> GSE
  // Applies rotations by the obliquity (ε) and the Sun’s ecliptic longitude (λ)
  // to reorient the inertial frame to one where the x‑axis points to the Sun and
  // the z‑axis is normal to the ecliptic.
  void t2(const std::array<int, 2>& itime, const std::array<double, 3>& x_in,
          std::array<double, 3>& x_out, int iverse) {
    double iyr, iday, ut;
    double tt0 = to_mjd(itime, iyr, iday, ut);
    double epsilon = (23.439 - 0.013 * tt0) * DEG2RAD;
    double m = (357.528 + 35999.05 * tt0 + 0.04107 * ut) * DEG2RAD;
    double cgamma = 280.46 + 36000.772 * tt0 + 0.04107 * ut;
    double lamdas = (cgamma + (1.915 - 0.0048 * tt0) * sin(m) + 0.02 * sin(2.0 * m)) * DEG2RAD;

    std::array<double, 3> temp;

    if (iverse == 1) {
      rotate_x(epsilon, x_in, temp);
      rotate_z(lamdas, temp, x_out);
    } else {
      rotate_z(double(iverse) * lamdas, x_in, temp);
      rotate_x(double(iverse) * epsilon, temp, x_out);
    }
  }

  // rotation matrix t3=<-psi,X>
  // GSE -> GSM
  // Rotates about the x‑axis by an angle (ψ) so that the magnetic dipole
  // projection falls into the desired plane.
  void t3(const std::array<int, 2>& itime, const std::array<double, 3>& x_in,
          std::array<double, 3>& x_out, int iverse) {
    double iyr, iday, ut;
    std::array<double, 3> q_c;

    get_q_c(itime, q_c);

    double psi = 0.0;

    if (q_c[2] == 0.0) {
      psi = -sign(1.5707963, q_c[1]);
    } else {
      psi = -atan(q_c[1] / q_c[2]);
    }
    rotate_x(double(iverse) * psi, x_in, x_out);
  }

  // transformation matrix t4=<-mu,Y>
  // GSM ->SM
  // Rotates about the y‑axis by an angle (μ) to align the magnetic dipole with
  // the z‑axis of the SM frame.
  void t4(const std::array<int, 2>& itime, const std::array<double, 3>& x_in,
          std::array<double, 3>& x_out, int iverse) {
    double iyr, iday, ut;
    std::array<double, 3> q_c;

    get_q_c(itime, q_c);

    double mu = -atan(q_c[0] / sqrt(q_c[1] * q_c[1] + q_c[2] * q_c[2]));
    rotate_y(double(iverse) * mu, x_in, x_out);
  }

  // transformation matrix t5=<phi-90,Y>*<lamda,Z>
  // GEO -> MAG
  // Rotates the geographic (Earth-fixed) coordinates by a geomagnetic longitude
  // (λ) and a latitude offset (φ – 90°) to produce a magnetic (geomagnetic)
  // coordinate system.
  void t5(const std::array<int, 2>& itime, const std::array<double, 3>& x_in,
          std::array<double, 3>& x_out, int iverse) {
    double pid2 = 3.1415927 / 2.0;

    double iyr = itime[0] / 1000;
    double iday = itime[0] - iyr * 1000;
    double ut = itime[1] / 1000.0 / 3600.0;
    double fracday = ut / 24.0;

    double rmjd = 45.0 + double(iyr - 1859) * 365.0 + (float((iyr - 1861) / 4) + 1.0) + double(iday)
                  - 1.0 + fracday;

    double factor = (rmjd - 46066.0) / 365.25;
    double phi = (78.8 + 4.283e-2 * factor) * DEG2RAD;
    double lamda = (289.1 - 1.413e-2 * factor) * DEG2RAD;

    std::array<double, 3> temp;
    if (iverse == 1) {
      rotate_z(lamda, x_in, temp);
      rotate_y(phi - pid2, temp, x_out);
    } else {
      rotate_y(iverse * (phi - pid2), x_in, temp);
      rotate_z(iverse * lamda, temp, x_out);
    }
  }

  void rotate_x(double angle, const std::array<double, 3>& x_in, std::array<double, 3>& x_out) {
    double c = cos(angle);
    double s = sin(angle);
    x_out[0] = x_in[0];
    x_out[1] = c * x_in[1] + s * x_in[2];
    x_out[2] = -s * x_in[1] + c * x_in[2];
  }

  void rotate_y(double angle, const std::array<double, 3>& x_in, std::array<double, 3>& x_out) {
    double c = cos(angle);
    double s = sin(angle);
    x_out[0] = c * x_in[0] + s * x_in[2];
    x_out[1] = x_in[1];
    x_out[2] = -s * x_in[0] + c * x_in[2];
  }

  void rotate_z(double angle, const std::array<double, 3>& x_in, std::array<double, 3>& x_out) {
    double c = cos(angle);
    double s = sin(angle);
    x_out[0] = c * x_in[0] + s * x_in[1];
    x_out[1] = -s * x_in[0] + c * x_in[1];
    x_out[2] = x_in[2];
  }

  // subroutine for computing the dipolar axis in GSE coordinates
  void get_q_c(const std::array<int, 2>& itime, std::array<double, 3>& q_c) {
    double iyr = itime[0] / 1000;
    double iday = itime[0] - iyr * 1000;
    double ut = itime[1] / 1000.0 / 3600.0;
    double fracday = ut / 24.0;

    double rmjd = 45.0 + double(iyr - 1859) * 365.0 + (float((iyr - 1861) / 4) + 1.0) + double(iday)
                  - 1.0 + fracday;

    double factor = (rmjd - 46066.0) / 365.25;
    double phi = (78.8 + 4.283e-2 * factor) * DEG2RAD;
    double lamda = (289.1 - 1.413e-2 * factor) * DEG2RAD;

    std::array<double, 3> temp, q_g;

    pol_to_cart(phi, lamda, 1.0, q_g);

    t1(itime, q_g, temp, -1);  // GEO -> GEI

    t2(itime, temp, q_c, 1);  // GEI -> GSE
  }

}  // namespace pecsim
