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

#include "lupnt/environment/occultation.h"

#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/constants.h"
#include "lupnt/data/kernels.h"
#include "lupnt/environment/body.h"

namespace lupnt {

  std::map<std::string, bool> Occultation::ComputeOccultation(
      Real epoch, const Vec3& r1, const Vec3& r2, Frame frame1, Frame frame2,
      const std::vector<BodyId>& bodies, const VecXd& height_atm, const VecXd& min_elevation,
      const bool use_elev_mask1 = false, const bool use_elev_mask2 = false) {
    // Check if the input vectors have the same size
    if (bodies.size() != height_atm.size())
      throw std::runtime_error("bodies and height_atm must have the same size");

    Vec3 r1_icrf = ConvertFrame(epoch, r1, frame1, Frame::ICRF);
    Vec3 r2_icrf = ConvertFrame(epoch, r2, frame2, Frame::ICRF);

    // Compute distance between line of sight and center of each body
    std::map<std::string, bool> vis;
    for (auto i = 0; i < bodies.size(); i++) {
      BodyData body_data = GetBodyData(bodies[i]);

      Vec3d rb = GetBodyPosVel(epoch, BodyId::SOLAR_SYSTEM_BARYCENTER, bodies[i], Frame::GCRF)
                     .head(3)
                     .cast<double>();
      Vec3d r12 = (r2_icrf - r1_icrf).cast<double>();  // 1->2
      Vec3d r1b = (rb - r1_icrf).cast<double>();       // 1->body
      Vec3d r2b = (rb - r2_icrf).cast<double>();       // 2->body

      // If the transmitter or receiver is inside the atmosphere, continue to elevation mask check
      // (Likely they are ground users)

      // Compute angle between (tx->Body center) and (tx->rx)
      double r1b_norm = r1b.norm();
      double r12_norm = r12.norm();
      double alpha_body = acos(r12.dot(r1b) / (r1b_norm * r12_norm));

      // Agent 1 or Agent 2 is ground user
      if ((r1b.norm() < height_atm(i)) || (r2b.norm() < height_atm(i))) {
        vis[body_data.name] = true;

        // Case 1: Agent 1 is ground user
        if ((use_elev_mask1) && (r1b.norm() < height_atm(i))) {
          // Angle between body->1->2
          if ((alpha_body - PI / 2.0) < min_elevation[i]) {
            vis[body_data.name] = false;
          }
        }
        // Case 2: Agent 2 is ground user
        else if ((use_elev_mask2) && (r2b.norm() < height_atm(i))) {
          // Angle between body->2->1
          // double b21 = acos(r2b.dot(-r12) / (r2b.norm() * r12.norm()));
          if ((alpha_body - PI / 2.0) < min_elevation[i]) {
            vis[body_data.name] = false;
          }
        }

        // Then, skip the rest of the occultation check for this body
        continue;
      }

      // Compute angle between (tx>Body center) and (tx->horizon)
      double R_body = GetBodyRadius(bodies[i]) + height_atm(i);
      double beta_body = asin(R_body / r1b_norm);

      // Compute occultation (alpha_body < beta and r12_norm > r1b_norm *
      // cos(beta))
      bool is_occult = (alpha_body < beta_body) && (r12_norm > r1b_norm * cos(beta_body));
      bool is_vis = !is_occult;

      vis[body_data.name] = is_vis;
    }

    vis["all"] = std::all_of(vis.begin(), vis.end(),
                             [](const std::pair<std::string, bool>& p) { return p.second; });
    return vis;
  }

  // MatX = func(real, Mat<-1, 3>, Mat<-1, 3>)
  std::vector<std::map<std::string, bool>> Occultation::ComputeOccultation(
      Real epoch, const Mat<-1, 3>& r1, const Mat<-1, 3>& r2, Frame frame1, Frame frame2,
      const std::vector<BodyId>& bodies, const VecXd& height_atm, const VecXd& min_elev,
      const bool use_elev_mask1 = false, const bool use_elev_mask2 = false) {
    if (!(r1.rows() == r2.rows() || r1.rows() == 1 || r2.rows() == 1))
      throw std::runtime_error(
          "r1 and r2 must have the same number of rows or one of them must have only one row");

    std::vector<std::map<std::string, bool>> vis;
    for (int i = 0; i < std::max(r1.rows(), r2.rows()); i++) {
      Vec3 r1_vec = r1.rows() == 1 ? r1.row(0) : r1.row(i);
      Vec3 r2_vec = r2.rows() == 1 ? r2.row(0) : r2.row(i);
      vis.push_back(ComputeOccultation(epoch, r1_vec, r2_vec, frame1, frame2, bodies, height_atm,
                                       min_elev, use_elev_mask1, use_elev_mask2));
    }
    return vis;
  }

  std::vector<std::map<std::string, bool>> Occultation::ComputeOccultation(
      const VecX& epoch, const Mat<-1, 3>& r1, const Mat<-1, 3>& r2, Frame frame1, Frame frame2,
      const std::vector<BodyId>& bodies, const VecXd& height_atm, const VecXd& min_elev,
      const bool use_elev_mask1 = false, const bool use_elev_mask2 = false) {
    if (!(epoch.size() == r1.rows() || r1.rows() == 1))
      throw std::runtime_error("epoch and r1 must have the same size or r1 must have only one row");
    if (!(epoch.size() == r2.rows() || r2.rows() == 1))
      throw std::runtime_error("epoch and r2 must have the same size or r2 must have only one row");

    std::vector<std::map<std::string, bool>> vis;
    for (int i = 0; i < epoch.size(); i++) {
      Vec3 r1_vec = r1.rows() == 1 ? r1.row(0) : r1.row(i);
      Vec3 r2_vec = r2.rows() == 1 ? r2.row(0) : r2.row(i);
      vis.push_back(ComputeOccultation(epoch(i), r1_vec, r2_vec, frame1, frame2, bodies, height_atm,
                                       min_elev, use_elev_mask1, use_elev_mask2));
    }
    return vis;
  }

}  // namespace lupnt
