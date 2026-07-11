// Satellite-to-Satellite Angles-Only Orbit Determination
// ------------------------------------------------------------------------------
// C++ counterpart of configs/sat_bearing_odts.yaml.
//
// Two `Spacecraft` agents share a `World` (the read-only physical environment):
//
//   * an OBSERVER spacecraft (low circular lunar orbit) hosts an `AnglesOdtsApp`,
//     which measures the unit line-of-sight (bearing) direction to the target and
//     runs an EKF that estimates its OWN orbit 6-state [r, v] from those bearings
//     (no clock -- angles are clock-independent);
//   * a TARGET spacecraft (elliptical lunar orbit) whose orbit/ephemeris is
//     treated as KNOWN (this is landmark-style navigation; the ill-posed
//     angles-only estimation of the target's RANGE is deliberately NOT attempted).
//
// The whole scenario lives in configs/sat_bearing_odts.yaml; the Simulation just
// holds the agents, the World, and an event queue, and Run() drains it (the app
// filters the whole arc in a single end-of-arc Monte-Carlo solve).

#include <yaml-cpp/yaml.h>

#include <iomanip>
#include <iostream>
#include <string>

#include "lupnt/applications/angles_odts/angles_odts_app.h"
#include "lupnt/lupnt.h"

using namespace lupnt;

int main(int argc, char** argv) {
  std::string config_path = (argc > 1) ? argv[1] : std::string("configs/sat_bearing_odts.yaml");

  YAML::Node config = YAML::LoadFile(config_path);
  Simulation sim(config);
  sim.Run();

  // The observer agent's application owns the angles-only OD solution.
  Agent* obs_agent = sim.GetAgent("observer");
  auto app = std::dynamic_pointer_cast<AnglesOdtsApp>(obs_agent->GetApplication());
  LUPNT_CHECK(app, "observer does not host an AnglesOdtsApp", "sat_bearing_odts");

  const VecXd& pos_rms = app->PositionErrorRms();
  const VecXd& vel_rms = app->VelocityErrorRms();
  const int N = static_cast<int>(pos_rms.size());

  std::cout << std::fixed << std::setprecision(3);
  std::cout << "Angles-only ODTS over " << N << " epochs, " << app->MonteCarloRuns()
            << " Monte-Carlo run(s)\n";
  std::cout << "Initial position error (RMS) : " << pos_rms(0) << " m\n";
  std::cout << "Final   position error (RMS) : " << app->FinalPositionErrorM() << " m\n";
  std::cout << "Tail    position error (RMS) : " << app->RmsPositionErrorM() << " m\n";
  std::cout << "Final   velocity error (RMS) : " << vel_rms(N - 1) * 1e3 << " mm/s\n";

  return 0;
}
