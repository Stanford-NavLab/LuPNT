// Example 7: Ground-Station Orbit Determination for a Lunar Satellite
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex7_groundstation_odts.ipynb.
//
// A lunar satellite in an Elliptical Lunar Frozen Orbit is tracked by the three
// 70 m Deep Space Network antennas (Goldstone, Canberra, Madrid) via two-way
// range and range-rate (Doppler). `GroundStationOdtsSimulation` propagates the
// truth trajectory, runs an elevation-mask visibility analysis, simulates noisy
// measurements over the visible passes, and recovers the orbit from a
// deliberately perturbed initial guess with an iterative batch (weighted
// least-squares) filter whose design matrix is built analytically from the
// autodiff state-transition matrix.

#include <iomanip>
#include <iostream>

#include "lupnt/lupnt.h"
#include "lupnt/simulations/ground_station_odts/ground_station_odts_simulation.h"

using namespace lupnt;

int main() {
  GroundStationOdtsConfig config;
  config.seed = 42;
  config.start_epoch_utc = "2026-01-01T00:00:00";
  config.duration_s = 36.0 * 3600.0;
  config.obs_interval_s = 300.0;
  config.elevation_mask_deg = 10.0;

  // Two-way range (10 m) and range-rate (1 mm/s) measurement noise.
  config.use_range = true;
  config.use_range_rate = true;
  config.range_sigma_m = 10.0;
  config.range_rate_sigma_mps = 1.0e-3;

  // A-priori error injected into the batch filter's starting guess.
  config.initial_position_sigma_m = 2000.0;
  config.initial_velocity_sigma_mps = 0.5;
  config.batch_use_analytic_jacobian = true;

  GroundStationOdtsSimulation sim(config);
  sim.Setup();
  sim.Precompute();
  sim.Run();
  const GroundStationOdtsResults& res = sim.GetResults();

  std::cout << "Ground stations: ";
  for (const std::string& n : res.station_names) std::cout << n << " ";
  std::cout << "\n";
  std::cout << "Simulated measurements: " << res.obs_epoch_index.size() << " over "
            << config.duration_s / 3600.0 << " h\n";
  std::cout << "Batch filter converged: " << (res.converged ? "yes" : "no") << " in "
            << res.num_iterations << " iterations\n\n";

  // Final solution error vs. truth initial state.
  const Vec3d dr = res.x0_estimated.head(3) - res.x0_true.head(3);
  const Vec3d dv = res.x0_estimated.tail(3) - res.x0_true.tail(3);
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "Initial-guess position error : "
            << (res.x0_initial_guess.head(3) - res.x0_true.head(3)).norm() << " m\n";
  std::cout << "Estimated  position error    : " << dr.norm() << " m\n";
  std::cout << "Estimated  velocity error    : " << dv.norm() * 1e3 << " mm/s\n";

  // Formal 1-sigma position uncertainty from the estimated covariance.
  const double sigma_pos = std::sqrt(res.covariance.topLeftCorner(3, 3).trace());
  std::cout << "Formal 1-sigma position (RSS): " << sigma_pos << " m\n";

  return 0;
}
