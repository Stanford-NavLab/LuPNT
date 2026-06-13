#include "lupnt/lupnt.h"
#include "src/filter_execution.h"
#include "src/gnssmeas_loader.h"
#include "src/simulation_config.h"
#include "src/state_meas_manager.h"

using namespace lupnt;
using std::cout;
using std::endl;
using namespace filtering_sim;

int main(int argc, char* argv[]) {
  // Simulation configuration
  SimulationConfig sim_config = LoadSimulationConfig("Plasmasphere_Delay_Filtering");

  // Open cache and output files
  H5Easy::File output_file = GetH5File(sim_config.output_path, true);

  // Copy config file to output directory if not exists
  if (!std::filesystem::exists(sim_config.output_dir / "main_config.yaml")) {
    std::filesystem::copy_file(sim_config.config_path, sim_config.output_dir / "main_config.yaml");
  }

  // Load Measurements from CSV files ----------------------------------------------------
  GnssMeasLoader gnss_meas_loader;
  gnss_meas_loader.set_min_cn0_dbhz(sim_config.min_cn0_dbhz);
  int file_i = 0;

  // First load receiver history
  gnss_meas_loader.load_rx_history(sim_config.rx_cache_path);

  for (const auto& csv_file : sim_config.meas_csv_files) {
    // Load and process each CSV file
    // Count Lines
    int sat_id = sim_config.meas_satids[file_i];
    std::string gnss = sim_config.meas_gnss[file_i];
    int signal = sim_config.meas_signals[file_i];
    double rz = sim_config.meas_rz[file_i];
    std::string orbit_dt_str = "norbit_" + std::to_string(static_cast<int>(sim_config.N_orbit))
                               + "_" + "dt_" + std::to_string(static_cast<int>(sim_config.dt))
                               + "s";

    std::string gnss_config_str = gnss + "_sat" + std::to_string(sat_id) + "_sig"
                                  + std::to_string(signal) + "_rz"
                                  + std::to_string(static_cast<int>(rz));
    if (!std::filesystem::exists(sim_config.output_dir / "gnss_cache" / orbit_dt_str
                                 / gnss_config_str)) {
      std::filesystem::create_directories(sim_config.output_dir / "gnss_cache" / orbit_dt_str
                                          / gnss_config_str);
    }

    std::cout << " " << std::endl;
    std::cout << "[Loading " << gnss << " L" << signal << " measurements]" << std::endl;
    std::string cache_file_name;
    std::string pco_clock_str = "ephemerr_" + std::to_string(int(sim_config.max_ephem_error_m));
    if (sim_config.correct_pco) {
      pco_clock_str += "_pco";
    }
    if (sim_config.correct_clock_bias) {
      pco_clock_str += "_clockbias";
    }
    cache_file_name
        = "cache_cn0_" + std::to_string(int(sim_config.min_cn0_dbhz)) + pco_clock_str + ".h5";

    std::filesystem::path cache_path
        = sim_config.output_dir / "gnss_cache" / orbit_dt_str / gnss_config_str / cache_file_name;
    gnss_meas_loader.load(csv_file, cache_path, sim_config.recompute_gnss, sim_config.correct_pco,
                          sim_config.correct_clock_bias, sim_config.max_ephem_error_m);
    file_i++;
  }

  gnss_meas_loader.construct_merged_posvel_histories();

  // Print statistics for delays
  gnss_meas_loader.print_delay_statistics();

  std::cout << "Finished loading all measurement CSV files." << std::endl;

  auto gnssmeas_data = gnss_meas_loader.get_data();
  cout << "Total measurements loaded: " << gnssmeas_data.size() << endl;
  VecXd tspan = gnss_meas_loader.get_tspan();
  MatX6d rx_posvel_mci = gnss_meas_loader.get_receiver_posvel_mci();
  cout << "Total unique time indices: " << tspan.size() << endl;

  // Run Monte Carlo Simulations --------------------------------------------------------
  auto pbar_all
      = Logger::GetProgressBar(sim_config.N_mc, "Running Monte-Carlo Simulations", "Main");

  bool verbose = true;

  for (int i = 0; i < sim_config.N_mc; i++) {
    std::cout << " " << std::endl;
    std::cout << std::string(200, '=') << std::endl;
    std::cout << "MC " << i + 1 << " / " << sim_config.N_mc << std::endl;
    std::cout << std::string(200, '=') << std::endl;
    std::cout << " " << std::endl;

    // First setup StateMeasManager
    StateMeasManager state_meas_manager
        = StateMeasManager(sim_config, gnss_meas_loader, i, sim_config.output_dir, verbose);
    // Run filtering
    filtering_sim::ExecuteFilter(state_meas_manager, i, verbose);
  }

  return 0;
}
