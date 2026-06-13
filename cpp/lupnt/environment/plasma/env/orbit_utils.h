/**
 * @file env/orbit_utils.h
 * @author Keidai Iiyama
 * @brief This file contains utility functions for orbit calculations.
 * @version 0.1
 * @date 2025-02-17
 */

#pragma once

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include "lupnt/environment/plasma/core/definitions.h"
#include "lupnt/environment/plasma/gcpm/constants_gcpm.h"

namespace pecsim {

  class Satellite {
  public:
    int id_;                     // Satellite ID (PRN)
    Vec6d posvel_;               // Position in ECEF coordinates [km]
    Vec6d coe_;                  // Classical orbital elements (a, e, i, Omega, w, M)
    double epoch_utc_;           // Epoch in UTC seconds (since J2000 epoch)
    double GM_ = 3.986004418e5;  // Gravitational constant for Earth [km^3/s^2]

    Satellite(int id, const Vec6d& coe, double epoch_utc, double GM = GM_EARTH);

    void propagate(double epoch_utc_new);
    Vec3d get_pos() const;           // Get position in ECI coordinates [km]
    Vec3d get_pos_epoch(double dt);  // Propagate position by time dt [s]
  };

  // Math utils
  double wrap2pi(double angle);  // Wrap angle to [0, 2*pi)
  Mat3d RotXd(double angle);
  Mat3d RotYd(double angle);  // Passive rotation matrix about the y-axis
  Mat3d RotZd(double angle);  // Passive rotation matrix about the z-axis

  // Anomalies
  double ecc2true(double E,
                  double e);  // Convert eccentric anomaly to true anomaly
  double ecc2mean(double E, double e);
  double mean2ecc(double M, double e);
  double mean2true(double M, double e);
  double true2ecc(double nu, double e);
  double true2mean(double nu, double e);
  Vec6d coe2cart(const Vec6d& coe, double GM);
  Vec6d cart2coe(const Vec6d& rv, double GM);

  // Orbit conversions
  Vec6d coe2cart(const Vec6d& coe, double GM);
  Vec6d cart2coe(const Vec6d& posvel, double GM);
  Vec6d propagate_coe(const Vec6d& coe, double GM,
                      double dt);  // Propagate COE

  // visibility
  bool compute_vis(const Vec3d& pos1, const Vec3d& pos2, double radius, double max_angle);
  double compute_min_altitude(const Vec3d& pos1, const Vec3d& pos2, double radius);
  Vec3d solve_lt(Satellite sat, const Vec3d& rx_pos,
                 double epoch_utc_rx);  // Solve light time for txpos

  // satellite
  std::vector<Satellite> setup_gnss_constellation(std::string filename);

}  // namespace pecsim
