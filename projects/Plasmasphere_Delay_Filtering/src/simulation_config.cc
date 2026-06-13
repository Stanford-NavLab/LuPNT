#include "simulation_config.h"

#include <lupnt/lupnt.h>

namespace filtering_sim {
  using namespace lupnt;

  SimulationConfig LoadSimulationConfig(const std::string& config_subdir) {
    SimulationConfig cfg;

    // Construct output path: output/project_name/v0
    const std::string project_name = "Plasmasphere_Delay_Filtering";
    const std::string version = "v0";
    std::filesystem::path base_output = GetBaseDir() / "output" / project_name / version;
    if (!std::filesystem::exists(base_output)) {
      std::filesystem::create_directories(base_output);
    }

    cfg.output_dir = base_output;
    cfg.cache_path = cfg.output_dir / "cache.h5";
    cfg.output_path = cfg.output_dir / "output.h5";
    cfg.config_path = GetBaseDir() / "projects" / config_subdir / "main_config.yaml";

    cfg.config = YAML::LoadFile(cfg.config_path);

    // Recompute flag
    cfg.recompute_gnss = cfg.config["recompute"]["gnss"].as<bool>(false);
    cfg.recompute_filter = cfg.config["recompute"]["filter"].as<bool>(false);

    std::cout << "[Recompute Settings] " << std::endl;
    std::cout << "  GNSS Measurement Cache: "
              << (cfg.recompute_gnss ? "Recompute" : "Use Existing Cache") << std::endl;
    std::cout << "  Filter Results: "
              << (cfg.recompute_filter ? "Recompute" : "Use Existing Results") << std::endl;
    std::cout << " " << std::endl;

    // Load time configuration
    auto time_cfg = cfg.config["time"];
    cfg.dt = time_cfg["dt"].as<double>();
    cfg.dt_prop = time_cfg["dt_prop"].as<double>();
    cfg.dt_raytrace = time_cfg["dt_raytrace"].as<double>();
    double N_orbit = time_cfg["N_orbit"].as<double>();
    cfg.N_orbit = N_orbit;

    // Load scenario configuration
    auto scenario_cfg = cfg.config["scenario"];
    cfg.case_name = scenario_cfg["case_name"].as<std::string>("test");

    // Monte Carlo simulations
    cfg.N_mc = scenario_cfg["N_mc"].as<int>(1);

    // CN0 Settings
    cfg.min_cn0_dbhz = scenario_cfg["min_cn0"].as<double>(18.0);

    // Nmber of smoothing iterations
    cfg.N_smooth_iter = cfg.config["filter"]["N_smooth_iter"].as<int>(0);

    // Correct Phase Center Offset (PCO)
    cfg.correct_pco = cfg.config["measurements"]["correct_pco"].as<bool>(false);
    cfg.correct_clock_bias = cfg.config["measurements"]["correct_clock_bias"].as<bool>(false);
    cfg.max_ephem_error_m = cfg.config["measurements"]["max_ephem_error_m"].as<double>(5.0);
    std::cout << "Correct for PCO in measurements: " << (cfg.correct_pco ? "Yes" : "No")
              << std::endl;
    std::cout << "Correct for Systematic Clock Bias in measurements: "
              << (cfg.correct_clock_bias ? "Yes" : "No") << std::endl;
    std::cout << "Max Ephem Error in measurements: " << cfg.max_ephem_error_m << " m" << std::endl;

    // Default time string
    cfg.t0_utc_string = "2025-03-01T12:00:00.0";

    // Set measurement csv file
    std::string norbit_dt_str = "norbit_" + std::to_string(static_cast<int>(N_orbit)) + "_dt_"
                                + std::to_string(static_cast<int>(cfg.dt)) + "s_dtrt_"
                                + std::to_string(static_cast<int>(cfg.dt_raytrace)) + "s";

    cfg.norbit_dt_str = norbit_dt_str;
    cfg.meas_csv_dir = GetOutputDir("iono_delay") / "raytrace_summary_pco" / norbit_dt_str;

    // Receiver Trajectory file
    cfg.sat_id = scenario_cfg["satellite_id"].as<int>(0);
    std::string rx_filename = "lcrns_sat" + std::to_string(cfg.sat_id) + "_norbit"
                              + std::to_string(static_cast<int>(N_orbit)) + "_dt"
                              + std::to_string(static_cast<int>(cfg.dt)) + ".h5";
    cfg.rx_cache_path = GetOutputDir("iono_delay") / "orbits" / "h5" / rx_filename;

    // Load all csv files in the directory
    std::cout << "Loading All Measurement CSV files from directory... " << std::endl;

    for (const auto& entry : std::filesystem::directory_iterator(cfg.meas_csv_dir)) {
      if (entry.path().extension() == ".csv") {
        // std::cout << "Found measurement CSV file: " << entry.path().string() << std::endl;

        // Get the filename
        std::string filename = entry.path().filename().string();
        // parse the filename (e.g.:
        // ionodata_raytrace_2025_03_01_12_00_00_sat_0_GALILEO_signal_5_rz12_50.0_kp_3.0.csv) and
        // get the sat number (int), gnss (string), signal (int), rz (double), kp (double) You can
        // store these values if needed for further processing
        std::istringstream iss(filename);
        std::string token;
        std::vector<std::string> tokens;
        while (std::getline(iss, token, '_')) {
          tokens.push_back(token);
        }
        // Filename:
        // ionodata_raytrace_2025_03_01_12_00_00_sat_0_GALILEO_signal_5_rz12_50.0_kp_3.0.csv
        // Indices:   0        1        2    3 4   5  6  7  8  9   10       11  12 13  14   15  16
        // Extract sat number, gnss, signal, rz
        int file_sat_id = std::stoi(tokens[9]);
        std::string gnss = tokens[10];
        int signal = std::stoi(tokens[12]);
        double rz = std::stod(tokens[14]);
        std::cout << "  Parsed sat_id: " << file_sat_id << ", gnss: " << gnss
                  << ", signal: " << signal << ", rz: " << rz << std::endl;

        // check if the extracted values are valid
        // gnss is in cfg.config["measurements"]["gnss"]
        // signal is in cfg.config["measurements"]["signal"]
        // rz is in cfg.config["measurements"]["rz"]
        auto valid_gnss = cfg.config["measurements"]["gnss"].as<std::vector<std::string>>();
        auto valid_signals = cfg.config["measurements"]["signal"].as<std::vector<int>>();
        auto valid_rz = cfg.config["measurements"]["rz"].as<std::vector<double>>();

        if (std::find(valid_gnss.begin(), valid_gnss.end(), gnss) == valid_gnss.end()) {
          std::cout << "Skipping file due to invalid GNSS: " << filename << std::endl;
          continue;
        }
        if (std::find(valid_signals.begin(), valid_signals.end(), signal) == valid_signals.end()) {
          std::cout << "Skipping file due to invalid signal: " << filename << std::endl;
          continue;
        }
        if (std::find(valid_rz.begin(), valid_rz.end(), rz) == valid_rz.end()) {
          std::cout << "Skipping file due to invalid Rz: " << filename << std::endl;
          continue;
        }

        // store the file path and parameters
        cfg.meas_csv_files.push_back(entry.path().string());
        cfg.meas_satids.push_back(file_sat_id);
        cfg.meas_gnss.push_back(gnss);
        cfg.meas_signals.push_back(signal);
        cfg.meas_rz.push_back(rz);
      }
    }

    // Logger::Info("Cache path: " + cfg.cache_path.string(), "Config");
    // Logger::Info("Output file: " + cfg.output_path.string(), "Config");
    // Logger::Info("Config file: " + cfg.config_path.string(), "Config");

    return cfg;
  }

}  // namespace filtering_sim
