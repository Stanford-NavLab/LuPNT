/**
 * @file raytrace_wrapper.cpp
 * @author Keidai Iiyama
 * @brief This file contains the wrapper functions for ray tracing
 * @version 0.1
 * @date 2025-02-17
 */

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <vector>

#if defined(_OPENMP) && defined(__has_include)
#  if __has_include(<omp.h>)
#    include <omp.h>
#    define PECSIM_HAS_OPENMP 1
#  endif
#elif defined(_OPENMP)
#  include <omp.h>
#  define PECSIM_HAS_OPENMP 1
#endif

#ifndef PECSIM_HAS_OPENMP
#  define omp_get_thread_num() 0
#  define omp_get_num_threads() 1
#endif

#include "lupnt/environment/plasma/tec/raytrace.h"
#include "lupnt/environment/plasma/tec/raytrace_wrapper.h"

namespace pecsim {

  RaytraceData load_raytrace_data(const std::string &filename) {
    RaytraceData data;

    std::ifstream file(filename);
    if (!file.is_open()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return data;
    }
    // Read header
    std::string header;
    std::getline(file, header);

    std::string line;
    std::regex vec_regex(R"(\[([^\]]+)\])");  // matches contents inside [...]

    while (std::getline(file, line)) {
      if (line.empty()) continue;
      std::stringstream ss(line);

      std::string field;
      std::vector<std::string> tokens;

      // Split by commas, but note that vectors contain commas only between
      // brackets
      while (std::getline(ss, field, ',')) {
        tokens.push_back(field);
      }

      if (tokens.size() < 12) {
        std::cerr << "Skipping malformed line: " << line << std::endl;
        continue;
      }

      int tidx = std::stoi(tokens[0]);
      data.tidxs.push_back(tidx);
      int row_full = std::stoi(tokens[1]);
      data.row_fulls.push_back(row_full);
      double epoch = std::stod(tokens[2]);
      data.epoch_utcs.push_back(epoch);
      double min_alt = std::stod(tokens[3]);
      data.min_alts.push_back(min_alt);

      auto parse_vec = [&](const std::string &s) -> Vec3d {
        std::smatch match;
        if (std::regex_search(s, match, vec_regex)) {
          std::stringstream vs(match[1].str());
          double a, b, c;
          vs >> a >> b >> c;
          return {a, b, c};
        }
        return {0.0, 0.0, 0.0};
      };

      Vec3d tx = parse_vec(tokens[4]);
      Vec3d rx = parse_vec(tokens[5]);
      data.tx_positions.push_back(tx);
      data.rx_positions.push_back(rx);

      data.total_delays.push_back(std::stod(tokens[6]));
      data.tecus.push_back(std::stod(tokens[7]));
      data.tec_delays.push_back(std::stod(tokens[8]));
      data.second_delays.push_back(std::stod(tokens[9]));
      data.third_delays.push_back(std::stod(tokens[10]));
      data.dist_bend_m.push_back(std::stod(tokens[11]));
      data.tec_delay_bend_m.push_back(std::stod(tokens[12]));
      data.max_sep_line_m.push_back(std::stod(tokens[13]));
      data.final_pos_err_m.push_back(std::stod(tokens[14]));
    }

    // std::cout << "Parsed " << data.epoch_utcs.size() << " rows." << std::endl;
    // if (!data.epoch_utcs.empty()) {
    //     std::cout << "Example:\n";
    //     std::cout << "  epoch_utc = " << data.epoch_utcs[0] << "\n";
    //     std::cout << "  tx = [" << data.tx_positions[0](0)<< ", "
    //               << data.tx_positions[0](1) << ", "
    //               << data.tx_positions[0](2) << "]\n";
    //     std::cout << "  total_delay = " << data.total_delays[0] << std::endl;
    // }

    return data;
  }

  void run_raytrace_batch(int job_id, int job_nums, const std::string epoch_str, const int sat_num,
                          const std::string gnss_const, const int freq_family, const double Rz12,
                          const double kp, bool correction, bool straight_ray,
                          const std::string datadir, const std::string savedir) {
    // --- One-time global IRI model setup (serial) ------------------------
    set_iri_model(IRIModel::IRI_2007);  // Set the IRI model to IRI-2007

    IRI2007Option op = IRI2007Option();
    op.R12 = Rz12;
    // If there is a function like set_iri_option(op), call it here.

    // Create a RayTraceConfig object
    RayTraceConfig config;
    if (freq_family == 1) {
      config.freq_Hz = freq_L1;
    } else if (freq_family == 2) {
      config.freq_Hz = freq_L2;
    } else if (freq_family == 5) {
      config.freq_Hz = freq_L5;
    } else {
      std::cerr << "Unsupported frequency family: " << freq_family << std::endl;
      return;
    }

    config.step_size = 20.0;
    config.correction = correction;
    config.fine_correction = false;
    config.cutoff_r = 4 * RE;
    config.gradn_dx = 1.0;
    config.integ_method = "RK4";
    config.kp = 3.0;
    config.correction_method = "neldermead";
    config.use_fortran_gcpm = false;
    config.corr_tol = 100.0;
    config.compute_higher_order = true;
    config.use_adaptive_step = true;
    config.straight_ray = straight_ray;

    bool debug_prop = false;
    bool debug_corr = false;

    // Load CSV file with ray path
    std::filesystem::path datadir_path(datadir);
    std::filesystem::path savedir_path(savedir);

    std::string sim_label_str_in = "sat_" + std::to_string(sat_num) + "_" + gnss_const + "_signal_"
                                   + std::to_string(freq_family);
    std::string sim_label_str_out = "sat_" + std::to_string(sat_num) + "_" + gnss_const + "_signal_"
                                    + std::to_string(freq_family) + "_rz12_"
                                    + std::to_string(int(Rz12)) + "_kp_" + std::to_string(int(kp));

    std::string filename
        = datadir_path / ("ionodata_short_" + epoch_str + "_" + sim_label_str_in + ".csv");
    std::cout << "Loading ray paths from: " << filename << std::endl;

    std::string datadir_name;
    if (straight_ray) {
      datadir_name = "raytrace_results_" + sim_label_str_out;
    } else {
      datadir_name = "raytrace_results_correct_" + sim_label_str_out;
    }

    std::string outdir = savedir_path / datadir_name;
    std::filesystem::create_directories(outdir);
    std::cout << "Output directory: " << outdir << std::endl;

    // Parse CSV file
    RaytraceData data = load_raytrace_data(filename);
    int num_rays = data.num_rays();

    if (num_rays == 0) {
      std::cerr << "No rays to process. Exiting." << std::endl;
      return;
    }

    std::cout << "Running raytrace for " << num_rays << " rays..." << std::endl;

    auto t_start = std::chrono::steady_clock::now();

    // Identify the idxs that requires computation for this job
    std::vector<int> needs_compute_idx;
    int total_needs_compute = 0;

    // Load Exisiting Data
    std::string outfilename;
    double min_alt_m, tecu, total_delay_m, dist_bend_m, tec_delay_m, tec_delay_bend_m,
        second_delay_m, third_delay_m, max_sep_line_m, final_pos_err_m;
    std::string outfile_headers
        = "total_delay_m,tecu,tec_delay_m,second_delay_m,third_delay_m,"
          "dist_bend_m,tec_delay_bend_m,max_sep_line_m,final_pos_err_m\n";

    for (int i = 0; i < num_rays; ++i) {
      outfilename = outdir + "/raytrace_result_" + std::to_string(i) + ".csv";

      if (straight_ray) {
        if ((data.total_delays[i] != 0.0) && (!std::filesystem::exists(outfilename))) {
          tecu = data.tecus[i];
          total_delay_m = data.total_delays[i];
          dist_bend_m = data.dist_bend_m[i];
          tec_delay_m = data.tec_delays[i];
          tec_delay_bend_m = data.tec_delay_bend_m[i];
          second_delay_m = data.second_delays[i];
          third_delay_m = data.third_delays[i];
          max_sep_line_m = data.max_sep_line_m[i];
          final_pos_err_m = data.final_pos_err_m[i];

          // Write to output file
          std::ofstream outfile(outfilename);
          outfile << outfile_headers;
          outfile << std::fixed << std::setprecision(6) << total_delay_m << "," << tecu << ","
                  << tec_delay_m << "," << second_delay_m << "," << third_delay_m << ","
                  << dist_bend_m << "," << tec_delay_bend_m << "," << max_sep_line_m << ","
                  << final_pos_err_m << "\n";
          outfile.close();
        } else if (!std::filesystem::exists(outfilename)) {
          needs_compute_idx.push_back(i);
          total_needs_compute++;
        }
      } else {  // not straight_ray
        if (!std::filesystem::exists(outfilename)) {
          if (data.min_alts[i] > 1500.0) {
            // when applying correction, limit the computation to rays with
            // min_alt < 1500 m to save time
            continue;
          }
          needs_compute_idx.push_back(i);
          total_needs_compute++;
        }
      }
    }

    // Devide the work into job_nums parts
    int start_idx = int(double(job_id * total_needs_compute) / double(job_nums));
    int end_idx = int(double(job_id + 1) * total_needs_compute / double(job_nums));
    if (job_id == job_nums - 1) {
      end_idx = total_needs_compute;
    }

    int num_sims = end_idx - start_idx;
    std::cout << "Total rays needing computation in this job: " << num_sims << std::endl;

    for (int i = start_idx; i < end_idx; ++i) {
      int idx = needs_compute_idx[i];
      double epoch_utc = data.epoch_utcs[idx];
      Vec3d pos_tx = data.tx_positions[idx];
      Vec3d pos_rx = data.rx_positions[idx];
      double total_delay_ref = data.total_delays[idx];
      min_alt_m = data.min_alts[idx];

      outfilename = outdir + "/raytrace_result_" + std::to_string(idx) + ".csv";

      // skip if output file already exists
      if (std::filesystem::exists(outfilename)) continue;

      if (!straight_ray) {
        if (min_alt_m < 1000.0) {
          // when applying correction, limit the computation to rays with
          // min_alt < 2000 m to save time
          config.correction = true;
        } else {
          config.correction = false;
        }
      }

      PathProfile pp = trace_ray(epoch_utc, pos_tx, pos_rx, config, debug_prop, debug_corr);

      tecu = pp.tecu;
      total_delay_m = pp.total_delay_m;
      dist_bend_m = pp.dist_bend_m;
      tec_delay_m = pp.tec_delay_m;
      tec_delay_bend_m = pp.tec_delay_bend_m;
      second_delay_m = pp.second_delay_m;
      third_delay_m = pp.third_delay_m;
      max_sep_line_m = pp.max_sep_line_m;
      final_pos_err_m = pp.corr_final_pos_err.norm() * 1000.0;

      std::ofstream outfile(outfilename);
      outfile << outfile_headers;
      outfile << std::fixed << std::setprecision(6) << total_delay_m << "," << tecu << ","
              << tec_delay_m << "," << second_delay_m << "," << third_delay_m << "," << dist_bend_m
              << "," << tec_delay_bend_m << "," << max_sep_line_m << "," << final_pos_err_m << "\n";
      outfile.close();

      // local progress index within this job
      int local_idx = i - start_idx;  // 0 .. num_sims-1
      int done = local_idx + 1;       // 1 .. num_sims

      int inv = 10;
      if (!straight_ray) {
        inv = 5;
      }

      if (done % inv == 0 || done == num_sims) {
        double pct = 100.0 * done / num_sims;
        auto t_now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::seconds>(t_now - t_start).count();
        double eta = (done > 0) ? elapsed * (num_sims - done) / done : 0.0;

        int elapsed_hr = int(elapsed / 3600);
        int elapsed_min = int((int(elapsed) % 3600) / 60);
        double elapsed_sec = int(elapsed) % 60;

        int eta_hr = int(eta / 3600);
        int eta_min = int((int(eta) % 3600) / 60);
        double eta_sec = int(eta) % 60;

        // single-line progress per job
        std::cout << "[Job " << job_id << "] Progress: " << done << " / " << num_sims << " ("
                  << std::fixed << std::setprecision(1) << pct << "%)"
                  << " | elapsed: " << elapsed_hr << "h " << elapsed_min << "m " << elapsed_sec
                  << "s"
                  << " | ETA: " << eta_hr << "h " << eta_min << "m " << eta_sec << "s   "
                  << std::endl;
      }
    }
    return;
  }

}  // namespace pecsim
