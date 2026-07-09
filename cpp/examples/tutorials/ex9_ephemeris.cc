// Example 9: Ephemeris and Almanac Fitting
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex9_ephemeris.ipynb.
//
// Studies the accuracy-vs-broadcast-datasize trade-off of two lunar-satellite
// orbit models fit to a numerically propagated truth trajectory:
//
//   * CartesianEphemeris -- a precise, short-validity Chebyshev model (the lunar
//     analogue of a GNSS broadcast ephemeris record).
//   * Almanac            -- a coarse, long-validity model (Algorithm 2 of
//     Iiyama & Gao) for whole-constellation acquisition.
//
// `EphemerisSimulation` fits both models over a sweep of fitting-window lengths,
// evaluating RMS / 95th-percentile RTN position and velocity error and the
// minimum broadcast bit budget needed to stay within a position tolerance.

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

int main() {
  EphemerisSimulationConfig config;
  config.start_epoch_utc = "2027-01-01T00:00:00";
  config.duration_days = 3.0;
  config.sample_dt_s = 60.0;

  // Truth force model: high-fidelity 20x20 lunar field + Earth/Sun third bodies.
  config.moon_gravity_degree = 20;
  config.moon_gravity_order = 20;
  config.include_earth = true;
  config.include_sun = true;

  // Fit and broadcast in the rotating Moon principal-axis frame (per Iiyama & Gao).
  config.output_frame = Frame::MOON_PA;

  // Sweep of fitting-window lengths and the required position accuracy.
  config.fit_window_minutes = {60.0, 120.0, 240.0, 480.0};
  config.num_windows = 6;
  config.datasize_precision_m = 0.01;

  EphemerisSimulation sim(config);
  sim.Setup();
  sim.Run();

  std::cout << "Ephemeris/Almanac fit study over " << config.duration_days << " days, "
            << "output frame MOON_PA\n";
  PrintTable("CartesianEphemeris (precise, short-validity):", sim.GetCartesianResults());
  PrintTable("Almanac (coarse, long-validity):", sim.GetAlmanacResults());

  return 0;
}
