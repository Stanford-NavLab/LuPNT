// Example 8: Distributed ISL + Lunar-Surface-Station ODTS
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex8_isl_odts.ipynb.
//
// A five-satellite lunar relay/navigation constellation (NASA LCRNS Reference
// Constellation 3.1), fully cross-linked, plus a small southern surface-station
// beacon network. The scenario is built from a single YAML file of PHYSICAL agents:
// each satellite is a `Spacecraft` (self-propagating orbit+clock truth) hosting a
// `SatelliteOdtsApp` onboard Schmidt-EKF that generates its own crosslinks + station
// aiding and exchanges consider-state with its neighbours; each `SurfaceStation`
// hosts a `StationBeaconSensor` feeding a `SurfaceStationManager`'s `GroundOdtsApp`
// (the centralized ground filter, no ISL). The shared environment (frame + truth
// force model) lives in the top-level `world:` block. `Simulation` runs the event
// loop; the per-epoch truth/estimate series are read off each app's accessors.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "lupnt/lupnt.h"

using namespace lupnt;

int main(int argc, char** argv) {
  const std::string config_path = argc > 1 ? argv[1] : "configs/isl_odts_distributed.yaml";

  Config cfg = YAML::LoadFile(config_path);
  Simulation sim(cfg);
  sim.Run();

  // Collect the physical satellite agents' onboard filters.
  std::vector<std::string> sat_names = {"SV-1", "SV-2", "SV-3", "SV-4", "SV-5"};
  std::vector<SatelliteOdtsApp*> onboard;
  for (const std::string& sn : sat_names) {
    auto* app = dynamic_cast<SatelliteOdtsApp*>(sim.GetAgent(sn)->GetApplication().get());
    LUPNT_CHECK(app, sn + " has no SatelliteOdtsApp", "ex8");
    onboard.push_back(app);
  }
  const int n_sat = static_cast<int>(onboard.size());
  const int N = static_cast<int>(onboard[0]->OwnEstimate().rows());
  std::cout << "Constellation: ";
  for (const std::string& n : sat_names) std::cout << n << " ";
  std::cout << "\n"
            << n_sat << " satellites, " << (n_sat - 1) << " crosslinks each, " << N << " epochs, "
            << n_sat << " onboard filters + 1 centralized ground filter\n\n";

  // Per-satellite onboard own position/clock error at the final epoch. OwnEstimate/TruthState
  // are the filter's own [r, v, clock_bias, clock_drift] 8-state ([N x 8]).
  std::cout << std::fixed << std::setprecision(3);
  for (int j = 0; j < n_sat; ++j) {
    const MatXd& est = onboard[j]->OwnEstimate();
    const MatXd& tru = onboard[j]->TruthState();
    const Vec3d r_true = tru.row(N - 1).head(3).transpose();
    const Vec3d r_est = est.row(N - 1).head(3).transpose();
    const double clk_err = tru(N - 1, 6) - est(N - 1, 6);
    std::cout << sat_names[j] << " onboard own pos err : " << (r_est - r_true).norm()
              << " m,  clock-bias err : " << clk_err * C << " m (range-equiv)\n";
  }

  // Centralized ground-filter (station pseudoranges only, no ISL) final position error.
  auto* ground = dynamic_cast<GroundOdtsApp*>(sim.GetAgent("gs_manager")->GetApplication().get());
  if (ground && ground->EstCentral().rows() == N) {
    std::cout << "\nCentralized ground filter (no ISL) final own pos err:\n";
    const MatXd& estc = ground->EstCentral();  // [N x 8*n_sat]
    for (int j = 0; j < n_sat; ++j) {
      const Vec3d r_true = onboard[j]->TruthState().row(N - 1).head(3).transpose();
      const Vec3d r_est = estc.row(N - 1).segment(8 * j, 3).transpose();
      std::cout << "  " << sat_names[j] << " : " << (r_est - r_true).norm() << " m\n";
    }
  }

  return 0;
}
