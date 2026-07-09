// Agent-based ground-station orbit determination.
//
// Loads a YAML scenario (an ELFO Satellite agent tracked by a GroundStation agent
// running a GroundStationOdtsApp), runs the Simulation event loop, and reports the
// batch orbit-determination result. This is the multi-agent counterpart of the
// monolithic GroundStationOdtsSimulation used by the ex7 notebook.
//
// Usage:
//   ./ex_ground_station_odts_app [scenario.yaml]
// Defaults to configs/ground_station_odts.yaml (run from the repository root).

#include <yaml-cpp/yaml.h>

#include <iostream>
#include <string>

#include "lupnt/applications/ground_station_odts_app.h"
#include "lupnt/lupnt.h"

using namespace lupnt;

int main(int argc, char** argv) {
  std::string config_path = (argc > 1) ? argv[1] : std::string("configs/ground_station_odts.yaml");
  std::cout << "Loading scenario: " << config_path << std::endl;

  YAML::Node config = YAML::LoadFile(config_path);
  Simulation sim(config);
  sim.Run();

  // Fetch the ground-station agent's application and report its ODTS solution.
  Agent* gs = sim.GetAgent("DSS14");
  auto app = std::dynamic_pointer_cast<GroundStationOdtsApp>(gs->GetApplication());
  LUPNT_CHECK(app, "DSS14 has no GroundStationOdtsApp", "ex_ground_station_odts_app");
  LUPNT_CHECK(app->HasSolved(), "GroundStationOdtsApp did not solve", "ex_ground_station_odts_app");

  Vec6d x0_true = app->X0True();
  Vec6d x0_est = app->X0Estimated();
  VecXd sigma = app->Covariance().diagonal().cwiseSqrt();
  double pos_err = (x0_est.head(3) - x0_true.head(3)).norm();
  double vel_err = (x0_est.tail(3) - x0_true.tail(3)).norm();

  std::cout << "\n=== Ground-station ODTS (agent-based) ===\n";
  std::cout << "Measurements collected: " << app->NumMeasurements() << "\n";
  std::cout << "Converged: " << (app->Converged() ? "yes" : "no") << "  (" << app->NumIterations()
            << " iterations)\n";
  std::cout << "Epoch state (MOON_CI):\n";
  std::cout << "  r0_true = " << x0_true.head(3).transpose() << " m\n";
  std::cout << "  r0_est  = " << x0_est.head(3).transpose() << " m\n";
  std::cout << "Final position error: " << pos_err << " m (formal 1-sigma " << sigma.head(3).norm()
            << " m)\n";
  std::cout << "Final velocity error: " << vel_err << " m/s (formal 1-sigma "
            << sigma.tail(3).norm() << " m/s)\n";

  return 0;
}
