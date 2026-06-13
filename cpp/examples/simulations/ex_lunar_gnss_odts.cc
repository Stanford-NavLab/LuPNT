// Command-line entry point for the Lunar GNSS ODTS scenario. The scenario can
// precompute links, run Monte Carlo filtering, or do both from the same config.
#include <filesystem>
#include <iostream>
#include <string>

#include "lupnt/lupnt.h"

namespace {

  void PrintUsage(const char* program) {
    std::cout << "Usage: " << program
              << " [--config PATH] [--precompute] [--run] [--precompute-and-run]\n";
  }

}  // namespace

int main(int argc, char** argv) {
  std::filesystem::path config_path = "projects/GNSS_Filtering/gnss_filtering_config.yaml";
  bool precompute = false;
  bool run = false;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--config") {
      if (i + 1 >= argc) {
        std::cerr << "--config requires a path\n";
        return 2;
      }
      config_path = argv[++i];
    } else if (arg == "--precompute") {
      precompute = true;
    } else if (arg == "--run") {
      run = true;
    } else if (arg == "--precompute-and-run") {
      precompute = true;
      run = true;
    } else if (arg == "-h" || arg == "--help") {
      PrintUsage(argv[0]);
      return 0;
    } else {
      std::cerr << "Unknown option: " << arg << "\n";
      PrintUsage(argv[0]);
      return 2;
    }
  }

  if (!precompute && !run) run = true;

  try {
    lupnt::LunarGnssODTSSimulation simulation(config_path);
    if (precompute) {
      simulation.Precompute();
      std::cout << "Lunar GNSS ODTS link precompute complete\n";
      std::cout << "  links: " << simulation.GetConfig().links_file << "\n";
    }
    if (run) {
      simulation.Run();
      const auto& summaries = simulation.GetSummaries();
      double mean_final_pos = 0.0;
      for (const auto& summary : summaries) mean_final_pos += summary.final_position_error_m;
      if (!summaries.empty()) mean_final_pos /= static_cast<double>(summaries.size());

      std::cout << "Lunar GNSS ODTS Monte Carlo complete\n";
      std::cout << "  runs: " << summaries.size() << "\n";
      std::cout << "  output: " << simulation.GetConfig().output_dir << "\n";
      std::cout << "  mean final position error [m]: " << mean_final_pos << "\n";
    }
  } catch (const std::exception& e) {
    std::cerr << "Lunar GNSS ODTS simulation failed: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
