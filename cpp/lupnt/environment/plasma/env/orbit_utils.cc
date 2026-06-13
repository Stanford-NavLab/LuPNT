/**
 * @file env/orbit_utils.h
 * @author Keidai Iiyama
 * @brief This file contains utility functions for orbit calculations.
 * @version 0.1
 * @date 2025-02-17
 */
#include "lupnt/environment/plasma/env/orbit_utils.h"

#include <cctype>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "lupnt/environment/plasma/core/string_file_utils.h"
#include "lupnt/environment/plasma/core/user_filepath.h"
#include "lupnt/environment/plasma/env/time_utils.h"
#include "lupnt/environment/plasma/gcpm/constants_gcpm.h"

namespace pecsim {

  double wrap2pi(double angle) {
    while (angle < 0.0) {
      angle += 2 * M_PI;
    }
    while (angle >= 2 * M_PI) {
      angle -= 2 * M_PI;
    }
    return angle;
  }

  /// @brief Passive rotation matrix about the x-axis
  /// @param angle Angle in radians
  /// @return Rotation matrix
  Mat3d RotXd(double angle) {
    double c = cos(angle);
    double s = sin(angle);
    Mat3d R{
        {1.0, 0.0, 0.0},
        {0.0, c, s},
        {0.0, -s, c},
    };
    return R;
  }

  /// @brief Passive rotation matrix about the y-axis
  /// @param angle Angle in radians
  /// @return Rotation matrix
  Mat3d RotYd(double angle) {
    double c = cos(angle);
    double s = sin(angle);
    Mat3d R{
        {c, 0.0, -s},
        {0.0, 1.0, 0.0},
        {s, 0.0, c},
    };
    return R;
  }

  /// @brief Passive rotation matrix about the z-axis
  /// @param angle Angle in radians
  /// @return Rotation matrix
  Mat3d RotZd(double angle) {
    double c = cos(angle);
    double s = sin(angle);
    Mat3d R{
        {c, s, 0.0},
        {-s, c, 0.0},
        {0.0, 0.0, 1.0},
    };
    return R;
  }

  /******************************************************************************
   * Conversion between different orbital anomalies
   * Eccentric anomaly (E), true anomaly (nu), and mean anomaly (M)
   * ***************************************************************************/

  double ecc2true(double E, double e) { return atan2(sqrt(1 - pow(e, 2)) * sin(E), cos(E) - e); }

  double ecc2mean(double E, double e) { return wrap2pi(E - e * sin(E)); }

  double mean2ecc(double M, double e) {
    double MM = wrap2pi(M);
    double EPS = 1e-10;  // Convergence criterion

    // Initial estimate of E
    double E = MM;
    double Eest = E - (E - e * sin(E) - MM) / (1.0 - e * cos(E));

    int max_itr = 100;
    int itr = 0;
    while ((abs(Eest - E) >= EPS) && (itr <= max_itr)) {
      E = Eest;
      Eest = E - (E - e * sin(E) - M) / (1.0 - e * cos(E));
      itr++;
    }
    E = Eest;

    return wrap2pi(E);
  }

  double mean2true(double M, double e) {
    double E = mean2ecc(M, e);
    double nu = ecc2true(E, e);
    return wrap2pi(nu);
  }

  double true2ecc(double nu, double e) {
    double E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(nu / 2));
    return E;
  }

  double true2mean(double nu, double e) {
    double E = true2ecc(nu, e);
    double M = ecc2mean(E, e);
    return M;
  }

  /******************************************************************************
   * Conversion between classical orbital elements (COE) and Cartesian coordinates
   * ***************************************************************************/

  Vec6d coe2cart(const Vec6d &coe, double GM) {
    double a = coe[0];      // Semi-major axis
    double e = coe[1];      // eccentricity
    double i = coe[2];      // Inclination
    double Omega = coe[3];  // Right ascension of ascending node
    double omega = coe[4];  // Argument of perigee
    double M = coe[5];      // mean anomaly

    double E = mean2ecc(M, e);  // eccentric anomaly
    double cosE = cos(E);
    double sinE = sin(E);

    double fac = sqrt((1.0 - e) * (1.0 + e));
    double R = a * (1.0 - e * cosE);  // Distance
    double V = sqrt(GM * a) / R;      // Velocity
    Vec3d r(a * (cosE - e), a * fac * sinE, 0.0);
    Vec3d v(-V * sinE, +V * fac * cosE, 0.0);

    Mat3d PQW = RotZd(-Omega) * RotXd(-i) * RotZd(-omega);
    r = PQW * r;
    v = PQW * v;
    Vec6d rv = Vec6d::Zero();
    rv.head<3>() = r;  // Position vector in ECEF coordinates
    rv.tail<3>() = v;  // Velocity vector in ECEF coordinates

    return rv;
  }

  Vec6d cart2coe(const Vec6d &rv, double GM) {
    Vec3d r = rv.head(3);  // Position
    Vec3d v = rv.tail(3);  // Velocity
    Vec3d h = r.cross(v);  // Areal velocity
    double H = h.norm();   // Angular momentum

    double Omega = atan2(h(0), -h(1));                        // Long. ascend. node
    double i = atan2(sqrt(h(0) * h(0) + h(1) * h(1)), h(2));  // Inclination
    double u = atan2(r(2) * H, -r(0) * h(1) + r(1) * h(0));   // Arg. of latitude
    double R = r.norm();                                      // Distance
    double a = 1.0 / (2.0 / R - v.squaredNorm() / GM);        // Semi-major axis

    double eCosE = 1.0 - R / a;              // e*cos(E)
    double eSinE = r.dot(v) / sqrt(GM * a);  // e*sin(E)
    double e2 = eCosE * eCosE + eSinE * eSinE;
    double e = sqrt(e2);             // eccentricity
    double E = atan2(eSinE, eCosE);  // eccentric anomaly

    double M = wrap2pi(E - eSinE);                          // mean anomaly
    double nu = atan2(sqrt(1.0 - e2) * eSinE, eCosE - e2);  // true anomaly
    double omega = wrap2pi(u - nu);                         // Arg. of perihelion
    Vec6d coe;
    coe << a, e, i, Omega, omega, M;  // Classical orbital elements
    return coe;
  }

  Vec6d propagate_coe(const Vec6d &coe, double GM, double dt) {
    // Propagate classical orbital elements by time dt
    double a = coe[0];      // Semi-major axis
    double e = coe[1];      // eccentricity
    double i = coe[2];      // Inclination
    double Omega = coe[3];  // Right ascension of ascending node
    double omega = coe[4];  // Argument of perigee
    double M = coe[5];      // mean anomaly

    double n = sqrt(GM / pow(a, 3));  // mean motion
    M = wrap2pi(M + n * dt);          // Update mean anomaly

    Vec6d new_coe;
    new_coe << a, e, i, Omega, omega, M;  // Updated classical orbital elements
    return new_coe;
  }

  /******************************************************************************
   * Satellite class implementation
   * This class represents a satellite in orbit, defined by its classical orbital
   * elements (COE) and provides methods to propagate its state and compute its
   * position at a given time.
   * ***************************************************************************/

  Satellite::Satellite(int id, const Vec6d &coe, double epoch_utc, double GM)
      : id_(id), coe_(coe), epoch_utc_(epoch_utc), GM_(GM) {
    // Convert COE to Cartesian position and velocity
    posvel_ = coe2cart(coe_, GM_);
  }

  void Satellite::propagate(double epoch_utc_new) {  // propagate to a new epoch in two-body
    // Propagate the satellite state to a new epoch
    Vec6d coe_new = propagate_coe(coe_, GM_, epoch_utc_new - epoch_utc_);

    // Update states
    coe_ = coe_new;
    posvel_ = coe2cart(coe_new, GM_EARTH);
    epoch_utc_ = epoch_utc_new;
  }

  Vec3d Satellite::get_pos_epoch(double epoch_utc_new) {
    // Propagate the satellite position by a time step dt
    Vec6d coe_new = propagate_coe(coe_, GM_EARTH, epoch_utc_new - epoch_utc_);
    Vec6d posvel = coe2cart(coe_new, GM_EARTH);
    return posvel.head<3>();  // Return new position
  }

  Vec3d Satellite::get_pos() const {
    // Get the position of the satellite in ECEF coordinates
    return posvel_.head<3>();
  }

  /******************************************************************************
   * Visbility function
   * ***************************************************************************/
  bool compute_vis(const Vec3d &tx_pos, const Vec3d &rx_pos, double radius, double max_angle) {
    Vec3d tx2rx = rx_pos - tx_pos;  // Vector from transmitter to receiver
    Vec3d tx2earth = -tx_pos;       // Vector from transmitter to Earth center
    double dot = tx2rx.dot(tx2earth);
    double mag_tx2rx = tx2rx.norm();
    double mag_tx2earth = tx2earth.norm();
    double angle = acos(dot / (mag_tx2rx * mag_tx2earth)) * 180 / M_PI;
    double angle_earth = asin(radius / mag_tx2earth) * 180 / M_PI;

    bool vis = true;
    if (angle > max_angle) {
      // std::cout << "Not within the mainlobe" << std::endl;
      vis = false;
    }
    if (angle < angle_earth) {
      // std::cout << "Blocked by Earth" << std::endl;
      vis = false;
    }

    return vis;
  }

  double compute_min_altitude(const Vec3d &pos1, const Vec3d &pos2, double radius) {
    Vec3d d = pos2 - pos1;
    double d_norm_sq = d.squaredNorm();

    if (d_norm_sq == 0.0) {
      // Degenerate case: same position
      return pos1.norm() - radius;
    }

    // Projection scalar t (clamped to [0,1] for line segment)
    double t = (-(pos1.dot(d))) / d_norm_sq;
    t = std::clamp(t, 0.0, 1.0);

    // Closest point on the segment
    Vec3d closest_point = pos1 + t * d;

    // Altitude = distance to Earth center - Earth radius
    return closest_point.norm() - radius;
  }

  Vec3d solve_lt(Satellite sat, const Vec3d &rx_pos, double epoch_utc) {
    // Solve light time for the transmission position
    Vec3d tx_pos = sat.get_pos_epoch(epoch_utc);
    Vec3d dir = (rx_pos - tx_pos).normalized();
    double L = (rx_pos - tx_pos).norm();
    double t_prop = L / C;             // Light time in seconds
    double t_tx = epoch_utc - t_prop;  // Transmission time in seconds since J2000

    // Recompute the transmission position at the transmission time
    tx_pos = sat.get_pos_epoch(t_tx);

    int n_iter = 2;  // Number of iterations for convergence
    for (int i = 0; i < n_iter; i++) {
      // Recompute the distance and light time
      L = (rx_pos - tx_pos).norm();
      t_prop = L / C;
      t_tx = epoch_utc - t_prop;

      // Update the transmission position
      tx_pos = sat.get_pos_epoch(t_tx);
    }

    return tx_pos;
  }

  /******************************************************************************
   * GNSS constellation setup
   * This function reads TLE data from a file and sets up the GNSS constellation
   * by creating Satellite objects for each satellite in the TLE data.
   * ***************************************************************************/
  std::vector<Satellite> setup_gnss_constellation(std::string filename) {
    int sat_num = 0;
    double epoch_utc = 0.0;
    std::vector<Satellite> satellites;

    for (auto tle : TLE::FromFile(filename)) {
      double sat_epoch = tle.epoch_utc;
      if (sat_epoch > epoch_utc) {
        epoch_utc = sat_epoch;  // Update epoch to the latest satellite epoch
      }

      double dt_epoch = sat_epoch - epoch_utc;
      double T = SECS_DAY / tle.mean_motion;

      // Classical orbital elements
      double a = pow((T * T * GM_EARTH) / (4.0 * PI * PI), 1.0 / 3.0);
      double e = tle.eccentricity;
      double inc = tle.inclination * DEG2RAD;
      double Omega = tle.raan * DEG2RAD;
      double w = tle.arg_perigee * DEG2RAD;
      double rad_per_sec = tle.mean_motion * 2 * PI / SECS_DAY;  // TLE mean motion is in revs/day
      double M = wrap2pi(tle.mean_anomaly * DEG2RAD + dt_epoch * rad_per_sec);  // [rad]

      // Convert to Cartesian
      Vec6d coe = {a, e, inc, Omega, w, M};

      Satellite sat = Satellite(tle.prn, coe,
                                sat_epoch);  // Create a Satellite object with PRN and posvel

      satellites.push_back(sat);
      sat_num++;
    }

    // propagate all satellites to the latest epoch
    for (auto &sat : satellites) {
      sat.propagate(epoch_utc);
    }

    std::cout << "Total number of satellites: " << sat_num << std::endl;
    std::cout << "Epoch UTC: " << epoch_utc << " seconds since J2000 epoch" << std::endl;

    return satellites;
  }

}  // namespace pecsim
