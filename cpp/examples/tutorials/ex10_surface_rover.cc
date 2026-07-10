// Example 10: Surface Rover Navigation -- IMU + LCRNS + DEM
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex10_surface_rover.ipynb.
//
// A lunar-surface rover near the south pole fuses a full strapdown IMU
// (accelerometer + gyroscope, Kalibr noise model), LCRNS/LANS pseudoranges from a
// five-satellite relay constellation, and a LOLA-DEM altitude constraint, with the
// IMU biases estimated online. The scenario is built from a single YAML file: a
// thin `Rover` agent hosts a `SurfaceRoverNavApp` (an error-state strapdown-INS
// EKF); the shared physical environment (frame, gravity, DEM site) lives in the
// top-level `world:` block. `pnt.Simulation` runs the event loop; the per-epoch
// truth/estimate/covariance series are read off the app.
//
// The last block re-runs the arc with the DEM constraint disabled to show how the
// terrain height ties down the weakly-observable vertical channel.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "lupnt/lupnt.h"

using namespace lupnt;

static double FinalUpError(const SurfaceRoverNavApp& app) {
  const MatXd& e = app.pos_err_enu();
  return std::abs(e(e.rows() - 1, 2));
}

int main(int argc, char** argv) {
  const std::string config_path = argc > 1 ? argv[1] : "configs/surface_rover_nav.yaml";

  // --- Run with the DEM altitude constraint enabled --------------------------
  Config cfg = YAML::LoadFile(config_path);
  Simulation sim(cfg);
  sim.Run();
  auto* app = dynamic_cast<SurfaceRoverNavApp*>(sim.GetAgent("Rover")->GetApplication().get());

  const int N = static_cast<int>(app->time_series().size());
  std::cout << "Site: " << app->site_id() << " (" << app->site_name() << ")\n";
  std::cout << "Ran " << N
            << " epochs, mean visible satellites: " << app->n_visible().cast<double>().mean()
            << "\n";
  std::cout << std::fixed << std::setprecision(2);
  std::cout << "Final 3D position error (DEM on) : " << app->pos_err_norm()(N - 1) << " m\n";

  // --- Ablation: same arc, DEM constraint disabled ---------------------------
  Config cfg_off = YAML::LoadFile(config_path);
  cfg_off["agents"]["Rover"]["application"]["enable_dem_constraint"] = false;
  Simulation sim_off(cfg_off);
  sim_off.Run();
  auto* app_off
      = dynamic_cast<SurfaceRoverNavApp*>(sim_off.GetAgent("Rover")->GetApplication().get());
  const int M = static_cast<int>(app_off->time_series().size());
  std::cout << "Final 3D position error (DEM off): " << app_off->pos_err_norm()(M - 1) << " m\n";
  std::cout << "Final Up-channel error   (DEM on) : " << FinalUpError(*app) << " m\n";
  std::cout << "Final Up-channel error   (DEM off): " << FinalUpError(*app_off) << " m\n";

  return 0;
}
