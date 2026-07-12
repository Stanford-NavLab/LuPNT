// Example 9: Ephemeris and LansAlmanac Fitting
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex9_ephemeris.ipynb.
//
// Studies the accuracy-vs-broadcast-datasize trade-off of two lunar-satellite
// orbit models fit to a numerically propagated truth trajectory:
//
//   * LansEphemeris -- a precise, short-validity Chebyshev model (the lunar
//     analogue of a GNSS broadcast ephemeris record).
//   * LansAlmanac            -- a coarse, long-validity model (Algorithm 2 of
//     Iiyama & Gao) for whole-constellation acquisition.
//
// The scenario is built from a single YAML file: a thin `EphemerisManager`
// coordinator agent hosts an `EphemerisApp`, and the shared environment (frame +
// a truth force model) lives in the top-level `world:` block. `Simulation` runs
// the event loop; the app's single scheduled step propagates the truth trajectory
// and, for each fitting-window length, fits both models over several windows,
// evaluating RMS / 95th-percentile RTN position and velocity error and the minimum
// broadcast bit budget needed to stay within a position tolerance. The per-window
// results are read off the app's `GetResults()`.

#include <iomanip>
#include <iostream>

#include "lupnt/lupnt.h"

using namespace lupnt;

static void PrintTable(const std::string& title, const std::vector<EphemerisWindowResult>& rows) {
  std::cout << "\n" << title << "\n";
  std::cout << "  window[min]  params  bits   pos_rms[m]  pos_p95[m]  vel_rms[mm/s]\n";
  std::cout << std::fixed;
  for (const EphemerisWindowResult& r : rows) {
    std::cout << "  " << std::setw(9) << std::setprecision(1) << r.fit_window_min << "  "
              << std::setw(6) << r.num_params << "  " << std::setw(5) << r.total_bits << "  "
              << std::setw(9) << std::setprecision(4) << r.pos_rms_m << "  " << std::setw(9)
              << r.pos_p95_m << "  " << std::setw(11) << std::setprecision(4) << r.vel_rms_mps * 1e3
              << "\n";
  }
}

int main(int argc, char** argv) {
  const std::string config_path = argc > 1 ? argv[1] : "configs/ephemeris.yaml";

  Config cfg = YAML::LoadFile(config_path);
  Simulation sim(cfg);
  sim.Run();
  auto* app = dynamic_cast<EphemerisApp*>(sim.GetAgent("EphemerisManager")->GetApplication().get());
  LUPNT_CHECK(app, "EphemerisManager agent has no EphemerisApp", "ex9");
  const EphemerisResults& res = app->GetResults();

  std::cout << "Ephemeris/LansAlmanac fit study over " << app->GetConfig().duration_days
            << " days, "
            << "output frame " << enum_name(app->GetConfig().output_frame) << "\n";
  PrintTable("LansEphemeris (precise, short-validity):", res.cartesian_results);
  PrintTable("LansAlmanac (coarse, long-validity):", res.almanac_results);

  return 0;
}
