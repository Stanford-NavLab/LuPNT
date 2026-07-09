// Example 11: Lunar Lander Navigation -- IMU + Altimeter + Crater Bearings + LunaNet
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex11_lander_navigation.ipynb.
//
// A lunar lander on powered descent to a south-pole site navigates with an
// error-state INS EKF that fuses a full IMU (Kalibr noise model), a nadir radar
// altimeter (height above DEM terrain), crater-landmark bearings (terrain-relative
// navigation against a synthetic crater map), and LunaNet/LANS pseudoranges, with
// the IMU biases estimated online. `RunLanderNav` loads the DEM, builds the crater
// map, generates the powered-descent truth trajectory, and drives the descent EKF.
//
// The final block re-runs with each aiding sensor disabled in turn to show what
// each one buys: craters -> horizontal, altimeter -> vertical, LunaNet -> absolute.

#include <iomanip>
#include <iostream>

#include "lupnt/lupnt.h"

using namespace lupnt;

// NASA LCRNS Reference Constellation 3.1 initial states, Frame::MOON_CI, [km].
static std::vector<LcrnsSatConfig> LunaNetConstellation() {
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

static double FinalPosError(LanderNavConfig cfg) {
  const LanderNavResults res = RunLanderNav(cfg);
  return res.pos_err_norm(static_cast<int>(res.pos_err_norm.size()) - 1);
}

int main() {
  LanderNavConfig config;
  config.seed = 42;
  config.start_epoch_utc = "2027-03-01T00:00:00";
  config.duration_s = 300.0;
  config.dt_s = 0.5;

  // Site01 "Connecting Ridge" near the lunar south pole.
  config.site_lat_deg = -89.45;
  config.site_lon_deg = 222.8;

  // Powered descent from 2 km altitude, 2.5 km downrange, to a hover at 15 m.
  config.descent_start_alt_m = 2000.0;
  config.descent_end_alt_m = 15.0;

  config.satellites = LunaNetConstellation();
  config.enable_altimeter = true;
  config.enable_craters = true;
  config.enable_lunanet = true;

  const LanderNavResults res = RunLanderNav(config);
  const int N = static_cast<int>(res.time_s.size());
  std::cout << "Site: " << res.site_id << " (" << res.site_name << ")\n";
  std::cout << "Ran " << N << " descent epochs\n";
  std::cout << std::fixed << std::setprecision(2);
  std::cout << "Touchdown altitude (truth) : " << res.alt_truth(N - 1) << " m\n";
  std::cout << "Final 3D position error    : " << res.pos_err_norm(N - 1) << " m\n";
  std::cout << "Final 3D velocity error    : " << res.vel_err_norm(N - 1) << " m/s\n\n";

  // --- Sensor ablation: disable one aiding source at a time ------------------
  std::cout << "Sensor ablation (final 3D position error):\n";
  std::cout << "  all sensors      : " << res.pos_err_norm(N - 1) << " m\n";
  {
    LanderNavConfig c = config;
    c.enable_craters = false;
    std::cout << "  no craters       : " << FinalPosError(c) << " m\n";
  }
  {
    LanderNavConfig c = config;
    c.enable_altimeter = false;
    std::cout << "  no altimeter     : " << FinalPosError(c) << " m\n";
  }
  {
    LanderNavConfig c = config;
    c.enable_lunanet = false;
    std::cout << "  no LunaNet        : " << FinalPosError(c) << " m\n";
  }

  return 0;
}
