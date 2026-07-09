// Example 10: Surface Rover Navigation -- IMU + LCRNS + DEM
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex10_surface_rover.ipynb.
//
// A lunar-surface rover near the south pole fuses a full strapdown IMU
// (accelerometer + gyroscope, Kalibr noise model), LCRNS/LANS pseudoranges from
// a five-satellite relay constellation, and a LOLA-DEM altitude constraint, with
// the IMU biases estimated online. `RunSurfaceNav` loads the DEM for the site
// (downloading from NASA PGDA if not cached), synthesizes the measurements, runs
// the strapdown-INS EKF, and returns the logged truth/estimate/covariance series.
//
// The last block re-runs the arc with the DEM constraint disabled to show how the
// terrain height ties down the weakly-observable vertical channel.

#include <iomanip>
#include <iostream>

#include "lupnt/lupnt.h"

using namespace lupnt;

// NASA LCRNS Reference Constellation 3.1 initial states, Frame::MOON_CI, [km].
static std::vector<LcrnsSatConfig> LcrnsConstellation() {
  struct Row {
    const char* name;
    Vec3d r_km, v_kmps;
  };
  const std::vector<Row> rows = {
      {"SV-1", {-198.931445, 386.023132, 3458.355412}, {-1.539867, 0.022243, -0.091059}},
      {"SV-2", {1089.980598, -2260.918271, -18984.562823}, {0.280301, -0.004049, 0.016575}},
      {"SV-3", {9187.665852, -260.470721, -8543.222759}, {0.273629, 0.417588, -0.313828}},
      {"SV-4", {6484.019002, 14329.883921, -7966.345594}, {-0.258433, 0.033951, 0.235264}},
      {"SV-5", {-5074.242314, 14473.022751, -8638.530033}, {-0.229931, -0.026926, -0.264381}},
  };
  std::vector<LcrnsSatConfig> sats;
  for (const Row& r : rows) {
    LcrnsSatConfig s;
    s.name = r.name;
    s.r0_m = r.r_km * 1e3;
    s.v0_mps = r.v_kmps * 1e3;
    sats.push_back(s);
  }
  return sats;
}

int main() {
  SurfaceNavConfig config;
  config.seed = 42;
  config.start_epoch_utc = "2027-03-01T00:00:00";
  config.duration_s = 1800.0;
  config.dt_s = 1.0;

  // Site01 "Connecting Ridge" near the lunar south pole.
  config.site_lat_deg = -89.45;
  config.site_lon_deg = 222.8;
  config.dem_half_width_m = 4000.0;

  config.satellites = LcrnsConstellation();
  config.pseudorange_sigma_m = 1.0;
  config.sise_m = 3.0;

  // --- Run with the DEM altitude constraint enabled --------------------------
  config.enable_dem_constraint = true;
  const SurfaceNavResults res = RunSurfaceNav(config);

  const int N = static_cast<int>(res.time_s.size());
  std::cout << "Site: " << res.site_id << " (" << res.site_name << ")\n";
  std::cout << "Ran " << N
            << " epochs, mean visible satellites: " << res.n_visible.cast<double>().mean() << "\n";
  std::cout << std::fixed << std::setprecision(2);
  std::cout << "Final 3D position error (DEM on) : " << res.pos_err_norm(N - 1) << " m\n";

  // --- Ablation: same arc, DEM constraint disabled ---------------------------
  config.enable_dem_constraint = false;
  const SurfaceNavResults res_no_dem = RunSurfaceNav(config);
  const int M = static_cast<int>(res_no_dem.time_s.size());
  std::cout << "Final 3D position error (DEM off): " << res_no_dem.pos_err_norm(M - 1) << " m\n";
  std::cout << "Final Up-channel error   (DEM on) : " << std::abs(res.pos_err_enu(N - 1, 2))
            << " m\n";
  std::cout << "Final Up-channel error   (DEM off): " << std::abs(res_no_dem.pos_err_enu(M - 1, 2))
            << " m\n";

  return 0;
}
