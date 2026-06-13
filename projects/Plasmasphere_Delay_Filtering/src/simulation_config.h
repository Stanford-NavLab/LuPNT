#pragma once

#include <lupnt/lupnt.h>
#include <yaml-cpp/yaml.h>

#include <filesystem>
#include <string>

namespace filtering_sim {
  using namespace lupnt;

  // Paths
  struct SimulationConfig {
    // Paths
    std::filesystem::path output_dir;
    std::filesystem::path cache_path;
    std::filesystem::path output_path;
    std::filesystem::path config_path;
    std::filesystem::path rx_cache_path;

    // User
    std::string case_name;  // Case name
    int sat_id;             // Satellite ID to simulate

    // Time configuration
    double dt;                  // [s] Simulation time step
    double dt_prop;             // [s] Propagation time step
    double dt_raytrace;         // [s] Raytracing time step
    double tf;                  // [s] Simulation final time
    double N_orbit;             // Number of orbits to simulate
    std::string t0_utc_string;  // Initial UTC time string

    std::string norbit_dt_str;  // String representing number of orbits and time steps

    // Measurements
    std::string meas_csv_dir;                 // Measurement CSV directory path
    std::vector<std::string> meas_csv_files;  // Measurement CSV file path
    std::vector<int> meas_satids;             // Satellite IDs corresponding to each CSV file
    std::vector<std::string> meas_gnss;       // GNSS constellations corresponding to each CSV file
    std::vector<int> meas_signals;            // Signal types corresponding to each CSV file
    std::vector<double> meas_rz;              // Rz values corresponding to each CSV file

    double min_cn0_dbhz;  // Minimum C/N0 for measurements

    int N_mc;  // Number of Monte Carlo simulations

    int N_smooth_iter;  // Number of smoothing iterations

    // Flags
    bool recompute_gnss = false;      // Whether to recompute GNSS measurement cache
    bool recompute_filter = false;    // Whether to recompute filter results
    bool correct_pco = false;         // Whether to correct for PCO in the measurements
    bool correct_clock_bias = false;  // Whether to correct for clock bias in the measurements
    double max_ephem_error_m
        = 5.0;  // Maximum ephemeris range error (in meters) to include measurements

    // YAML config node
    YAML::Node config;
  };

  SimulationConfig LoadSimulationConfig(const std::string& config_subdir
                                        = "2025_ION_Nav_SurfaceStation");
}  // namespace filtering_sim
