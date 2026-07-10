// Example 8: Distributed ISL + Lunar-Surface-Station ODTS
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex8_isl_odts.ipynb.
//
// A five-satellite lunar relay/navigation constellation (NASA LCRNS Reference
// Constellation 3.1), fully cross-linked, plus a small southern surface-station
// beacon network. The scenario is built from a single YAML file: a thin
// `IslOdtsManager` coordinator agent hosts an `IslOdtsCoordinatorApp`, and the
// shared environment (frame + truth force model) lives in the top-level `world:`
// block. `Simulation` runs the event loop; each scheduled epoch the coordinator app
// propagates every satellite's truth trajectory, simulates two-way range/Doppler +
// time/frequency-transfer crosslinks between every pair, runs one onboard Schmidt-EKF
// per satellite in parallel (each estimating its own [r, v, clock_bias, clock_drift]
// while carrying the others as consider states), serves a one-way surface-station
// pseudorange/Doppler to visible satellites, runs a centralized ground filter, and
// periodically exchanges consider-state between the parallel filters. The per-epoch
// truth/estimate/covariance series are read off the coordinator app's `GetResults()`.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "lupnt/lupnt.h"

using namespace lupnt;

int main(int argc, char** argv) {
  const std::string config_path = argc > 1 ? argv[1] : "configs/isl_odts.yaml";

  Config cfg = YAML::LoadFile(config_path);
  Simulation sim(cfg);
  sim.Run();
  auto* app
      = dynamic_cast<IslOdtsCoordinatorApp*>(sim.GetAgent("IslManager")->GetApplication().get());
  LUPNT_CHECK(app, "IslManager agent has no IslOdtsCoordinatorApp", "ex8");
  const IslOdtsResults& res = app->GetResults();

  const int N = static_cast<int>(res.t_s.size());
  const int n_sat = static_cast<int>(res.satellite_names.size());
  std::cout << "Constellation: ";
  for (const std::string& n : res.satellite_names) std::cout << n << " ";
  std::cout << "\n"
            << n_sat << " satellites, " << (n_sat - 1) << " crosslinks each, " << N << " epochs, "
            << n_sat << " parallel filters\n\n";

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

  // Surface-station beacon coverage: satellite-epochs with at least one station in view.
  int visible_count = 0;
  for (int k = 0; k < N; ++k)
    for (int j = 0; j < n_sat; ++j)
      if (res.station_visible(k, j) > 0.0) visible_count++;
  std::cout << "\nStation beacon coverage: " << visible_count << " / " << (N * n_sat)
            << " satellite-epochs in view\n";

  // Centralized ground-filter final position error (station pseudoranges only, no ISL).
  if (res.est_central.size() > 0) {
    std::cout << "Centralized ground filter (no ISL) final own pos err:\n";
    for (int j = 0; j < n_sat; ++j) {
      const Vec3d r_true = res.truth_states[j].row(N - 1).head(3).transpose();
      const Vec3d r_est = res.est_central.row(N - 1).segment(8 * j, 3).transpose();
      std::cout << "  " << res.satellite_names[j] << " : " << (r_est - r_true).norm() << " m\n";
    }
  }

  return 0;
}
