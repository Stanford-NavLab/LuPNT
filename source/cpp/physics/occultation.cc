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

#include "lupnt/physics/occultation.h"

#include "lupnt/core/constants.h"
#include "lupnt/data/kernels.h"
#include "lupnt/physics/body.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {

  std::map<std::string, bool> Occultation::ComputeOccultationGnss(const Vec3d tx_eci,
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
  std::map<std::string, bool> Occultation::ComputeOccultation(Real epoch, const Vec3& r1,
                                                              const Vec3& r2, Frame cs1, Frame cs2,
                                                              const std::vector<NaifId>& bodies,
                                                              const VecXd& atm_h) {
    // Check if the input vectors have the same size
    assert(bodies.size() == atm_h.size() && "bodies and atm_h must have the same size");

    // Convert to ICRF
    Vec3 r1_icrf = ConvertFrame(epoch, r1, cs1, Frame::ICRF);
    Vec3 r2_icrf = ConvertFrame(epoch, r2, cs2, Frame::ICRF);

    // Compute distance between line of sight and center of each body
    std::map<std::string, bool> vis;

    for (size_t i = 0; i < bodies.size(); i++) {
      BodyData bodydata = GetBodyData(bodies[i]);

      // Vec6 GetBodyPosVel(const real t_tai, NaifId center, NaifId target);
      Vec3d rb = GetBodyPosVel(epoch, NaifId::SOLAR_SYSTEM_BARYCENTER, bodies[i], Frame::ICRF)
                     .head(3)
                     .cast<double>();
      Vec3d r12 = (r2_icrf - r1_icrf).cast<double>();  // 1->2
      Vec3d r1b = (rb - r1_icrf).cast<double>();       // 1->body
      Vec3d r2b = (rb - r2_icrf).cast<double>();       // 2->body

      // If the transmitter or receiver is inside the atmosphere, ignore
      // ocuultation computation (Likely they are ground users)
      if ((r1b.norm() < atm_h(i)) || (r2b.norm() < atm_h(i))) {
        vis[bodydata.name] = true;
        continue;
      }

      // Compute angle between (tx->Body center) and (tx->rx)
      double r1b_norm = r1b.norm();
      double r12_norm = r12.norm();

      double alpha_body = acos(r12.dot(r1b) / (r1b_norm * r12_norm));

      // Compute angle between (tx>Body center) and (tx->horizon)
      double R_body = GetBodyRadius(bodies[i]) + atm_h(i);
      double beta_body = asin(R_body / r1b_norm);

      // Compute occultation (alpha_body < beta and r12_norm > r1b_norm *
      // cos(beta))
      bool is_vis = (alpha_body < beta_body) && (r12_norm > r1b_norm * cos(beta_body));

      vis[bodydata.name] = is_vis;
    }

    vis["all"] = std::all_of(vis.begin(), vis.end(),
                             [](const std::pair<std::string, bool>& p) { return p.second; });
    return vis;
  }

  // MatX = func(real, Mat<-1, 3>, Mat<-1, 3>)
  std::vector<std::map<std::string, bool>> Occultation::ComputeOccultation(
      Real epoch, const Mat<-1, 3>& r1, const Mat<-1, 3>& r2, Frame cs1, Frame cs2,
      const std::vector<NaifId>& bodies, const VecXd& atm_h) {
    assert((r1.rows() == r2.rows() || r1.rows() == 1 || r2.rows() == 1) &&
         "r1 and r2 must have the same number of rows or one of them must have "
         "only one row");
    std::vector<std::map<std::string, bool>> vis;

    for (int i = 0; i < std::max(r1.rows(), r2.rows()); i++) {
      Vec3 r1_vec = r1.rows() == 1 ? r1.row(0) : r1.row(i);
      Vec3 r2_vec = r2.rows() == 1 ? r2.row(0) : r2.row(i);
      vis.push_back(ComputeOccultation(epoch, r1_vec, r2_vec, cs1, cs2, bodies, atm_h));
    }
    return vis;
  }

  std::vector<std::map<std::string, bool>> Occultation::ComputeOccultation(
      const VecX& epoch, const Mat<-1, 3>& r1, const Mat<-1, 3>& r2, Frame cs1, Frame cs2,
      const std::vector<NaifId>& bodies, const VecXd& atm_h) {
    assert((epoch.size() == r1.rows() || r1.rows() == 1)
           && "epoch and r1 must have the same size or r1 must have only one row");
    assert((epoch.size() == r2.rows() || r2.rows() == 1)
           && "epoch and r2 must have the same size or r2 must have only one row");
    std::vector<std::map<std::string, bool>> vis;
    for (int i = 0; i < epoch.size(); i++) {
      Vec3 r1_vec = r1.rows() == 1 ? r1.row(0) : r1.row(i);
      Vec3 r2_vec = r2.rows() == 1 ? r2.row(0) : r2.row(i);
      vis.push_back(ComputeOccultation(epoch(i), r1_vec, r2_vec, cs1, cs2, bodies, atm_h));
    }
    return vis;
  }

}  // namespace lupnt
