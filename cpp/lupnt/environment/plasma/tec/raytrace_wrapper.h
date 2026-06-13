/**
 * @file raytrace_wrapper.h
 * @author Keidai Iiyama
 * @brief This file contains the wrapper functions for ray tracing
 * @version 0.1
 * @date 2025-02-17
 *
 */

#pragma once

#include <string>
#include <vector>

#include "lupnt/environment/plasma/core/definitions.h"

namespace pecsim {

  struct RaytraceData {
    std::vector<int> tidxs;                // Transmitter indices
    std::vector<int> row_fulls;            // Full row indices
    std::vector<double> epoch_utcs;        // UTC time of each epoch (s)
    std::vector<double> min_alts;          // Minimum altitudes along the ray paths (m)
    std::vector<Vec3d> tx_positions;       // Transmitter positions [m]
    std::vector<Vec3d> rx_positions;       // Receiver positions [m]
    std::vector<double> total_delays;      // Total signal delay [m]
    std::vector<double> tecus;             // Total electron content [TECU]
    std::vector<double> tec_delays;        // TEC-induced delay [m]
    std::vector<double> second_delays;     // Second-order ionospheric delay [m]
    std::vector<double> third_delays;      // Third-order ionospheric delay [m]
    std::vector<double> dist_bend_m;       // Geometric path bending distance [m]
    std::vector<double> tec_delay_bend_m;  // TEC delay due to bending [m]
    std::vector<double> max_sep_line_m;    // Max separation between straight and bent paths [m]
    std::vector<double> final_pos_err_m;   // Final position error [m]
    int num_rays() const { return epoch_utcs.size(); }
  };

  RaytraceData load_raytrace_data(const std::string& filename);

  void run_raytrace_batch(int job_id, int job_nums, const std::string epoch_str, const int sat_idx,
                          const std::string gnss_const, const int freq_family, const double Rz12,
                          const double kp, bool correction, bool straight_ray,
                          const std::string datadir, const std::string savedir);

}  // namespace pecsim
