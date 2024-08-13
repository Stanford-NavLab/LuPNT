/**
 * @file occultation.h
 * @author Stanford NAV LAB
 * @brief Signal blockage model
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/measurements/occultation.h"

#include "lupnt/core/constants.h"
#include "lupnt/data/kernels.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {

  std::map<std::string, bool> Occultation::ComputeOccultation(const Vec3d tx_eci,
                                                              const Vec3d tx_mci,
                                                              const Vec3d rx_eci,
                                                              const Vec3d rx_mci,
                                                              const std::string tx_planet) {
    Vec3d tx2usr = rx_eci - tx_eci;
    double tx2usr_norm = tx2usr.norm();

    // COMPUTE EARTH OCCULTATION

    // Compute angle between (tx->Earth center) and (tx->rx)
    Vec3d tx2earth = -tx_eci;
    double tx2earth_norm = tx2earth.norm();
    double alpha_earth = acos(tx2usr.dot(tx2earth) / (tx2earth_norm * tx2usr_norm));

    bool occ_earth, occ_atmos, occ_ionos;
    if (tx_planet != "EARTH") {
      // Compute angle between (tx>Earth center) and (tx->horizon)
      double beta_earth = asin(R_EARTH / tx2earth_norm);
      double beta_atmos = asin(r_atmos_ / tx2earth_norm);
      double beta_ionos = asin(r_ionos_ / tx2earth_norm);

      // Compute occultation (alpha_earth < beta and tx2usr_norm >
      // tx2hor_norm)
      occ_earth = ((alpha_earth < beta_earth) && (tx2usr_norm > tx2earth_norm * cos(beta_earth)));
      occ_atmos = ((alpha_earth < beta_atmos) && (tx2usr_norm > tx2earth_norm * cos(beta_atmos)));
      occ_ionos = ((alpha_earth < beta_ionos) && (tx2usr_norm > tx2earth_norm * cos(beta_ionos)));
    } else {
      occ_earth = (alpha_earth < (PI / 2.0 + min_elevation_));
      occ_atmos = true;
      occ_ionos = true;
    }

    // COMPUTE MOON OCCULTATION

    // Compute angle between (tx->Moon center) and (tx->rx)
    Vec3d tx2moon = -tx_mci;
    double tx2moon_norm = tx2moon.norm();
    double alpha_moon = acos(tx2moon.dot(tx2usr) / (tx2moon_norm * tx2usr_norm));

    bool occ_moon;
    if (tx_planet != "MOON") {
      // Compute angle between (tx->Moon center) and (tx->horizon)
      double beta_moon = asin(R_MOON / tx2moon_norm);

      // Compute occultation (alpha_moon < beta and tx2usr_norm > tx2hor_norm)
      occ_moon = ((alpha_moon < beta_moon) && (tx2usr_norm > tx2moon_norm * cos(beta_moon)));
    } else {
      occ_moon = (alpha_moon < (PI / 2.0 + min_elevation_));
    }

    return {{"earth", occ_earth}, {"atmos", occ_atmos}, {"ionos", occ_ionos}, {"moon", occ_moon}};
  }

  // Return distance from the line of sight to the center of each body
  // VecX = func(real, Vec3, Vec3)
  VecX Occultation::ComputeOccultation(Real epoch, const Vec3& r1, const Vec3& r2, Frame cs1,
                                       Frame cs2, const std::vector<NaifId>& bodies) {
    // Convert to ICRF
    Vec3 r1_icrf = ConvertFrame(epoch, r1, cs1, Frame::ICRF);
    Vec3 r2_icrf = ConvertFrame(epoch, r2, cs2, Frame::ICRF);

    // Compute distance between line of sight and center of each body
    VecX distances(bodies.size());
    for (size_t i = 0; i < bodies.size(); i++) {
      // Vec6 GetBodyPosVel(const real t_tai, NaifId center, NaifId target);
      Vec3 rb = GetBodyPosVel(epoch, NaifId::SOLAR_SYSTEM_BARYCENTER, bodies[i]).head(3);
      Vec3 r12 = r2_icrf - r1_icrf;  // 1->2
      Vec3 r1b = rb - r1_icrf;       // 1->body
      Vec3 r2b = rb - r2_icrf;       // 2->body
      if (r1b.norm() <= r2b.norm()) {
        if (r1b.dot(r12) <= 0) {
          // Body is behind 1
          distances(i) = r1b.norm();
        } else {
          // Compute distance to line of sight
          Vec3 r12_unit = r12.normalized();
          Vec3 rblos = r1b.dot(r12_unit) * r12_unit - r1b;  // body->los
          distances(i) = rblos.norm();
        }
      } else {
        if (r2b.dot(-r12) <= 0) {
          // Body is behind 2
          distances(i) = r2b.norm();
        } else {
          // Compute distance to line of sight
          Vec3 r12_unit = -r12.normalized();
          Vec3 rblos = r2b.dot(r12_unit) * r12_unit - r2b;  // body->los
          distances(i) = rblos.norm();
        }
      }
    }
    return distances;
  }

  // MatX = func(real, Mat<-1, 3>, Mat<-1, 3>)
  MatX Occultation::ComputeOccultation(Real epoch, const Mat<-1, 3>& r1, const Mat<-1, 3>& r2,
                                       Frame cs1, Frame cs2, const std::vector<NaifId>& bodies) {
    assert((r1.rows() == r2.rows() || r1.rows() == 1 || r2.rows() == 1) &&
    "r1 and r2 must have the same number of rows or one of them must have "
    "only one row");
    MatX distances(std::max(r1.rows(), r2.rows()), bodies.size());
    for (int i = 0; i < std::max(r1.rows(), r2.rows()); i++) {
      Vec3 r1_vec = r1.rows() == 1 ? r1.row(0) : r1.row(i);
      Vec3 r2_vec = r2.rows() == 1 ? r2.row(0) : r2.row(i);
      distances.row(i) = ComputeOccultation(epoch, r1_vec, r2_vec, cs1, cs2, bodies);
    }
    return distances;
  }

  // MatX = func(VecX, Mat<-1, 3>, Mat<-1, 3>)
  MatX Occultation::ComputeOccultation(const VecX& epoch, const Mat<-1, 3>& r1,
                                       const Mat<-1, 3>& r2, Frame cs1, Frame cs2,
                                       const std::vector<NaifId>& bodies) {
    assert((epoch.size() == r1.rows() || r1.rows() == 1)
           && "epoch and r1 must have the same size or r1 must have only one row");
    assert((epoch.size() == r2.rows() || r2.rows() == 1)
           && "epoch and r2 must have the same size or r2 must have only one row");
    MatX distances(epoch.size(), bodies.size());
    for (int i = 0; i < epoch.size(); i++) {
      Vec3 r1_vec = r1.rows() == 1 ? r1.row(0) : r1.row(i);
      Vec3 r2_vec = r2.rows() == 1 ? r2.row(0) : r2.row(i);
      distances.row(i) = ComputeOccultation(epoch(i), r1_vec, r2_vec, cs1, cs2, bodies);
    }
    return distances;
  }

}  // namespace lupnt
