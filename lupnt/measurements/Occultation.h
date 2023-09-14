/**
 * @file Occultation.h
 * @author Guillem Casadesus Vila
 * @brief Occultation between Earth and Moon
 * @version 0.1
 * @date 2023-03-07
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <lupnt/core/Constants.h>

#include <autodiff/forward/real/eigen.hpp>
#include <cmath>
#include <map>

namespace LPT {

class Occultation {
 private:
  static constexpr double min_elevation_ = 10.0 * RAD_PER_DEG;
  static constexpr double r_atmos_ = R_EARTH + 100.0;  // atmospheric mask
  static constexpr double r_ionos_ = R_EARTH + 965.0;  // ionospheric mask

 public:
  // TODO: Generalize this
  Occultation();
  ~Occultation();

  /**
   * @brief Compute occultation between a tx and a rx
   *
   * @param tx_eci   The tx position in the Earth-centered inertial
   * @param tx_mci   The tx position in the Moon-centered inertial
   * @param rx_eci      The rx position in the Earth-centered inertial frame
   * @param rx_mci      The rx position in the Moon-centered inertial frame
   * @param tx_planet    The planet surface the transmitter is on
   * @return std::map<std::string, bool>
   */
  static std::map<std::string, bool> ComputeOccultation(
      const Eigen::Vector3d tx_eci, const Eigen::Vector3d tx_mci,
      const Eigen::Vector3d rx_eci, const Eigen::Vector3d rx_mci,
      const std::string tx_planet) {
    Eigen::Vector3d tx2usr = rx_eci - tx_eci;
    double tx2usr_norm = tx2usr.norm();

    // COMPUTE EARTH OCCULTATION

    // Compute angle between (tx->Earth center) and (tx->rx)
    Eigen::Vector3d tx2earth = -tx_eci;
    double tx2earth_norm = tx2earth.norm();
    double alpha_earth =
        acos(tx2usr.dot(tx2earth) / (tx2earth_norm * tx2usr_norm));

    bool occ_earth, occ_atmos, occ_ionos;
    if (tx_planet != "EARTH") {
      // Compute angle between (tx>Earth center) and (tx->horizon)
      double beta_earth = asin(R_EARTH / tx2earth_norm);
      double beta_atmos = asin(r_atmos_ / tx2earth_norm);
      double beta_ionos = asin(r_ionos_ / tx2earth_norm);

      // Compute occultation (alpha_earth < beta and tx2usr_norm >
      // tx2hor_norm)
      occ_earth = ((alpha_earth < beta_earth) &&
                   (tx2usr_norm > tx2earth_norm * cos(beta_earth)));
      occ_atmos = ((alpha_earth < beta_atmos) &&
                   (tx2usr_norm > tx2earth_norm * cos(beta_atmos)));
      occ_ionos = ((alpha_earth < beta_ionos) &&
                   (tx2usr_norm > tx2earth_norm * cos(beta_ionos)));
    } else {
      occ_earth = (alpha_earth < (M_PI / 2.0 + min_elevation_));
      occ_atmos = true;
      occ_ionos = true;
    }

    // COMPUTE MOON OCCULTATION

    // Compute angle between (tx->Moon center) and (tx->rx)
    Eigen::Vector3d tx2moon = -tx_mci;
    double tx2moon_norm = tx2moon.norm();
    double alpha_moon =
        acos(tx2moon.dot(tx2usr) / (tx2moon_norm * tx2usr_norm));

    bool occ_moon;
    if (tx_planet != "MOON") {
      // Compute angle between (tx->Moon center) and (tx->horizon)
      double beta_moon = asin(R_MOON / tx2moon_norm);

      // Compute occultation (alpha_moon < beta and tx2usr_norm > tx2hor_norm)
      occ_moon = ((alpha_moon < beta_moon) &&
                  (tx2usr_norm > tx2moon_norm * cos(beta_moon)));
    } else {
      occ_moon = (alpha_moon < (M_PI / 2.0 + min_elevation_));
    }

    return {{"earth", occ_earth},
            {"atmos", occ_atmos},
            {"ionos", occ_ionos},
            {"moon", occ_moon}};
  }
};
}  // namespace LPT
