// Example 11: Lunar Lander Navigation -- IMU + Altimeter + Crater Bearings + LunaNet
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex11_lander_navigation.ipynb.
//
// A lunar lander on powered descent to a south-pole site navigates with an
// error-state INS MEKF that fuses a full IMU (Kalibr noise model), a nadir radar
// altimeter (height above DEM terrain), crater-landmark bearings (terrain-relative
// navigation against a synthetic crater map), and LunaNet/LANS pseudoranges, with
// the IMU biases estimated online. The scenario is built from a single YAML file:
// a thin `Lander` agent hosts a `LanderNavApp`; the shared environment (frame,
// gravity, DEM site) lives in the top-level `world:` block. `pnt.Simulation` runs
// the descent; results are read off the app.
//
// The final block re-runs with each aiding sensor disabled in turn to show what
// each one buys: craters -> horizontal, altimeter -> vertical, LunaNet -> absolute.

#include <iomanip>
#include <iostream>
#include <string>

#include "lupnt/lupnt.h"

using namespace lupnt;

static double FinalPosError(Config cfg) {
  Simulation sim(cfg);
  sim.Run();
  auto* app = dynamic_cast<LanderNavApp*>(
      sim.GetAgent("Lander")->GetApplicationByName("LanderNavApp").get());
  const VecXd& e = app->pos_err_norm();
  return e(e.size() - 1);
}

int main(int argc, char** argv) {
  const std::string config_path = argc > 1 ? argv[1] : "configs/lander_nav.yaml";

  Config cfg = YAML::LoadFile(config_path);
  Simulation sim(cfg);
  sim.Run();
  auto* app = dynamic_cast<LanderNavApp*>(
      sim.GetAgent("Lander")->GetApplicationByName("LanderNavApp").get());

  const int N = static_cast<int>(app->time_series().size());
  std::cout << "Site: " << app->site_id() << " (" << app->site_name() << ")\n";
  std::cout << "Ran " << N << " descent epochs\n";
  std::cout << std::fixed << std::setprecision(2);
  std::cout << "Touchdown altitude (truth) : " << app->alt_truth()(N - 1) << " m\n";
  std::cout << "Final 3D position error    : " << app->pos_err_norm()(N - 1) << " m\n";
  std::cout << "Final 3D velocity error    : " << app->vel_err_norm()(N - 1) << " m/s\n\n";

  // --- Sensor ablation: disable one aiding source at a time ------------------
  std::cout << "Sensor ablation (final 3D position error):\n";
  std::cout << "  all sensors      : " << app->pos_err_norm()(N - 1) << " m\n";
  {
    Config c = YAML::LoadFile(config_path);
    c["agents"]["Lander"]["application"]["enable_craters"] = false;
    std::cout << "  no craters       : " << FinalPosError(c) << " m\n";
  }
  {
    Config c = YAML::LoadFile(config_path);
    c["agents"]["Lander"]["application"]["enable_altimeter"] = false;
    std::cout << "  no altimeter     : " << FinalPosError(c) << " m\n";
  }
  {
    Config c = YAML::LoadFile(config_path);
    c["agents"]["Lander"]["application"]["enable_lunanet"] = false;
    std::cout << "  no LunaNet        : " << FinalPosError(c) << " m\n";
  }

  return 0;
}
