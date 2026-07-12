// Example 6: Lunar GNSS Sidelobe Orbit Determination & Timing (ODTS)
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex6_gnss_odts.ipynb.
//
// A lunar-orbiting receiver navigates from Earth-GNSS sidelobe signals: sidelobe
// pseudorange plus time-differenced carrier phase (TDCP) drive an EKF that
// estimates the receiver orbit, clock (bias + drift), and SRP coefficient, with
// a GCPM/IRI plasmaspheric-delay stage on the truth signal. As with the ex7/ex8
// tutorials, the whole scenario is assembled from a single YAML file: a physical
// `receiver` Spacecraft agent hosts a `LunarGnssOdtsApp` whose scheduled step
// runs the entire Monte-Carlo ODTS body (bit-identical to the former
// `LunarGnssODTSSimulation`). `Simulation::Run()` drives it; the per-seed
// summaries (final / RMS position, velocity, and clock errors) are read off the
// app.
//
// Data note: a real run needs the SP3 precise ephemerides
// (data/LuPNT_data/ephemeris/gnsslibpy/sp3) and, ideally, the Stage 1 link +
// Stage 2 plasma-delay precompute caches (see python/examples/ex6_precompute.py
// and docs/pages/sp3_download). Absent the data the program prints an
// informative note and exits cleanly, so it always compiles and runs; with the
// data present, run it as `./ex6_gnss_odts` from the repository root.

#include <yaml-cpp/yaml.h>

#include <iomanip>
#include <iostream>
#include <string>

#include "lupnt/applications/lunar_gnss_odts/lunar_gnss_odts_app.h"
#include "lupnt/lupnt.h"

using namespace lupnt;

int main(int argc, char** argv) {
  const std::string config_path = argc > 1 ? argv[1] : std::string("configs/lunar_gnss_odts.yaml");

  try {
    Config cfg = YAML::LoadFile(config_path);
    Simulation sim(cfg);
    std::cout << "Running the lunar GNSS ODTS EKF (LunarGnssOdtsApp) ...\n";
    sim.Run();

    // The receiver Spacecraft hosts the ODTS coordinator application.
    auto app
        = std::dynamic_pointer_cast<LunarGnssOdtsApp>(sim.GetAgent("receiver")->GetApplication());
    LUPNT_CHECK(app, "receiver agent has no LunarGnssOdtsApp", "ex6_gnss_odts");

    const std::vector<LunarGnssODTSSummary>& summaries = app->GetSummaries();
    if (summaries.empty()) {
      std::cout << "Run completed but produced no summaries.\n";
      return 0;
    }

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nMonte-Carlo seeds: " << summaries.size() << "\n";
    const LunarGnssODTSSummary& s = summaries.front();
    std::cout << "Epochs processed         : " << s.num_epochs << "\n";
    std::cout << "Final position error     : " << std::setw(9) << s.final_position_error_m
              << " m\n";
    std::cout << "Final velocity error     : " << std::setw(9) << s.final_velocity_error_mps * 1e3
              << " mm/s\n";
    std::cout << "Final clock bias error   : " << std::setw(9) << s.final_clock_bias_error_m * 1e2
              << " cm\n";
    std::cout << "Final clock drift error  : " << std::setw(9)
              << s.final_clock_drift_error_mps * 1e3 << " mm/s\n";
    std::cout << "RMS position error       : " << std::setw(9) << s.rms_position_error_m << " m\n";
    std::cout << "RMS velocity error       : " << std::setw(9) << s.rms_velocity_error_mps * 1e3
              << " mm/s\n";

    if (summaries.size() > 1) {
      double pos = 0.0, rms = 0.0;
      for (const LunarGnssODTSSummary& e : summaries) {
        pos += e.final_position_error_m;
        rms += e.rms_position_error_m;
      }
      const double n = static_cast<double>(summaries.size());
      std::cout << "\nSeed-averaged final position error: " << pos / n << " m\n";
      std::cout << "Seed-averaged RMS   position error: " << rms / n << " m\n";
    }

  } catch (const std::exception& e) {
    std::cout << "\n[skipped] ODTS inputs unavailable: " << e.what() << "\n"
              << "This example needs the SP3 precise ephemerides under\n"
              << "data/LuPNT_data/ephemeris/gnsslibpy/sp3 and (ideally) the Stage 1/2\n"
              << "precompute caches. Build them with python/examples/ex6_precompute.py;\n"
              << "see docs/pages/sp3_download. The program compiled and ran; it skipped\n"
              << "the ODTS run because the inputs are absent.\n";
    return 0;
  }

  return 0;
}
