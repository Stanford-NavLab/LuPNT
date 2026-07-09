// Example 8: Distributed ISL + Lunar-Surface-Station ODTS
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex8_isl_odts.ipynb.
//
// A five-satellite lunar relay/navigation constellation (NASA LCRNS Reference
// Constellation 3.1), fully cross-linked. `IslOdtsSimulation` propagates every
// satellite's truth trajectory, simulates two-way range/Doppler crosslink
// measurements between every pair, and runs one onboard Schmidt-EKF per
// satellite in parallel: each estimates its own [r, v, clock_bias, clock_drift]
// while carrying the others as consider states. A lunar surface station serves
// the satellites one at a time (round-robin) with a one-way pseudorange, adding
// absolute position and clock observability. Every consider_exchange_interval_s
// the filters swap their own estimate (mean + covariance) to refresh each
// other's consider blocks.

#include <cmath>
#include <iomanip>
#include <iostream>

#include "lupnt/lupnt.h"

using namespace lupnt;

int main() {
  // --- NASA LCRNS Reference Constellation 3.1 (Ryden & Volle, GSFC), Table 2 --
  // Cartesian initial states in the Moon-centred ICRF (== Frame::MOON_CI), [km].
  // SV-1 is the hub. Clock offsets are representative (not in the reference doc).
  struct LcrnsRow {
    const char* name;
    Vec3d r_km, v_kmps;
    double clock_bias_s, clock_drift_sps;
  };
  const std::vector<LcrnsRow> lcrns = {
      {"SV-1", {-198.931445, 386.023132, 3458.355412}, {-1.539867, 0.022243, -0.091059}, 0.0, 0.0},
      {"SV-2",
       {1089.980598, -2260.918271, -18984.562823},
       {0.280301, -0.004049, 0.016575},
       10e-6,
       1e-10},
      {"SV-3",
       {9187.665852, -260.470721, -8543.222759},
       {0.273629, 0.417588, -0.313828},
       -6e-6,
       -8e-11},
      {"SV-4",
       {6484.019002, 14329.883921, -7966.345594},
       {-0.258433, 0.033951, 0.235264},
       4e-6,
       5e-11},
      {"SV-5",
       {-5074.242314, 14473.022751, -8638.530033},
       {-0.229931, -0.026926, -0.264381},
       -8e-6,
       -3e-11},
  };

  IslOdtsConfig config;
  config.seed = 42;
  config.start_epoch_utc = "2027-03-01T00:00:00";
  config.duration_s = 6.0 * 3600.0;
  config.dt_s = 60.0;
  for (const LcrnsRow& s : lcrns) {
    IslOdtsSatelliteConfig sat;
    sat.name = s.name;
    sat.r0_m = s.r_km * 1e3;
    sat.v0_mps = s.v_kmps * 1e3;
    sat.clock_bias_s = s.clock_bias_s;
    sat.clock_drift_sps = s.clock_drift_sps;
    config.satellites.push_back(sat);
  }

  config.moon_gravity_degree_truth = 16;
  config.moon_gravity_order_truth = 16;

  // Two-way crosslink measurement noise.
  config.range_sigma_m = 1.0;
  config.range_rate_sigma_mps = 1.0e-3;

  // Each filter knows its own state to ~200 m but has only a coarse ~500 m broadcast
  // prior of its neighbors -- the gap the inter-agent exchange closes.
  config.initial_position_sigma_m = 200.0;
  config.consider_position_sigma_m = 500.0;
  config.consider_velocity_sigma_mps = 0.05;
  config.consider_clock_bias_sigma_s = 1.0e-6;
  config.consider_clock_drift_sigma_sps = 1.0e-9;

  // Rotating lunar surface station (near the south pole) + 10-min filter exchange.
  config.surface_station.enabled = true;
  config.surface_station.latitude_deg = -89.9;
  config.surface_station.pseudorange_sigma_m = 5.0;
  config.consider_exchange_interval_s = 600.0;

  IslOdtsSimulation sim(config);
  sim.Setup();
  sim.Run();
  const IslOdtsResults& res = sim.GetResults();

  const int N = static_cast<int>(res.t_s.size());
  const int n_sat = static_cast<int>(res.satellite_names.size());
  std::cout << "Constellation: ";
  for (const std::string& n : res.satellite_names) std::cout << n << " ";
  std::cout << "\n" << n_sat << " satellites, " << (n_sat - 1) << " crosslinks each, " << N
            << " epochs, " << n_sat << " parallel filters\n\n";

  // Per-satellite own position/clock error at the final epoch. est[j][:, 8*j:8*j+8]
  // is filter j's own state (global sat order).
  std::cout << std::fixed << std::setprecision(3);
  for (int j = 0; j < n_sat; ++j) {
    const Vec3d r_true = res.truth_states[j].row(N - 1).head(3).transpose();
    const Vec3d r_est = res.est[j].row(N - 1).segment(8 * j, 3).transpose();
    const double clk_err = res.truth_states[j](N - 1, 6) - res.est[j](N - 1, 8 * j + 6);
    std::cout << res.satellite_names[j] << " final own pos err : " << (r_est - r_true).norm()
              << " m,  clock-bias err : " << clk_err * C << " m (range-equiv)\n";
  }

  // Surface-station rotation summary.
  int served_count = 0;
  for (int k = 0; k < N; ++k)
    if (res.served_sat_idx(k) >= 0) served_count++;
  std::cout << "\nStation served " << served_count << " / " << (N - 1) << " epochs\n";

  return 0;
}
