// Example 7: Ground-Station Orbit Determination for a Lunar Satellite
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex7_groundstation_odts.ipynb.
//
// A lunar satellite in an Elliptical Lunar Frozen Orbit is tracked by the three
// 70 m Deep Space Network antennas (Goldstone, Canberra, Madrid) via two-way
// range and range-rate (Doppler). The scenario is assembled from first-class
// agents that share a `World` (the read-only physical environment):
//
//   * a `Satellite` agent (`sat`) self-propagates the truth trajectory using the
//     World's force model;
//   * three `GroundStation` agents each run a `GroundStationTrackingApp` (a
//     sensor) that generates its own visibility-gated range/range-rate
//     observations and pushes them to the manager;
//   * a `GroundStationManager` agent (`gs_manager`) runs the centralized
//     `GroundStationManagerApp`, which aggregates all stations' observations and
//     recovers the orbit with an analytic batch (weighted least-squares) filter
//     plus a square-root information filter / smoother.
//
// The whole scenario lives in configs/ground_station_odts.yaml; the Simulation
// just holds the agents, the World, and an event queue, and Run() drains it.

#include <yaml-cpp/yaml.h>

#include <iomanip>
#include <iostream>
#include <string>

#include "lupnt/applications/ground_station/ground_station_manager_app.h"
#include "lupnt/lupnt.h"

using namespace lupnt;

int main(int argc, char** argv) {
  std::string config_path = (argc > 1) ? argv[1] : std::string("configs/ground_station_odts.yaml");

  YAML::Node config = YAML::LoadFile(config_path);
  Simulation sim(config);
  sim.Run();

  // The manager agent's application owns the centralized OD solution.
  Agent* mgr_agent = sim.GetAgent("gs_manager");
  auto mgr = std::dynamic_pointer_cast<GroundStationManagerApp>(mgr_agent->GetApplication());
  LUPNT_CHECK(mgr && mgr->HasSolved(), "GroundStationManagerApp did not solve",
              "ex7_groundstation_odts");

  std::cout << "Ground stations: ";
  for (const std::string& n : mgr->StationNames()) std::cout << n << " ";
  std::cout << "\n";
  std::cout << "Aggregated measurements: " << mgr->NumMeasurements() << "\n";
  std::cout << "Batch filter converged: " << (mgr->Converged() ? "yes" : "no") << " in "
            << mgr->NumIterations() << " iterations\n\n";

  const Vec6d x0_true = mgr->X0True();
  const Vec6d x0_est = mgr->X0Estimated();
  const Vec3d dr = x0_est.head(3) - x0_true.head(3);
  const Vec3d dv = x0_est.tail(3) - x0_true.tail(3);
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "Initial-guess position error : "
            << (mgr->X0InitialGuess().head(3) - x0_true.head(3)).norm() << " m\n";
  std::cout << "Estimated  position error    : " << dr.norm() << " m\n";
  std::cout << "Estimated  velocity error    : " << dv.norm() * 1e3 << " mm/s\n";

  const double sigma_pos = std::sqrt(mgr->Covariance().topLeftCorner(3, 3).trace());
  std::cout << "Formal 1-sigma position (RSS): " << sigma_pos << " m\n";

  return 0;
}
