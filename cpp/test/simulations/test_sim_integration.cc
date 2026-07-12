// Integration tests: build a `lupnt::Simulation` from a repo YAML config, shorten
// its duration in-memory (never touching the config files on disk), run the event
// loop, and assert basic sanity off the hosted applications.
//
// These are the highest-leverage coverage tests in the suite: a single Run()
// exercises World + Simulation + Agent + Application + measurements + filters +
// dynamics together. Every config here is self-contained (analytic orbits +
// synthetic measurements) and needs NO Earthdata/SP3 downloads. Durations are cut
// to a few epochs so each test finishes in ~seconds and is deterministic (the
// scenarios carry fixed seeds).
//
// The config path is resolved off `LUPNT_DATA_PATH` (== <root>/data/LuPNT_data),
// so it is independent of the test binary's working directory.

#include <lupnt/applications/ephemeris/ephemeris_app.h>
#include <lupnt/applications/ground_station/ground_station_manager_app.h>
#include <lupnt/applications/lander/lander_nav_app.h>
#include <lupnt/applications/lunar_sat_odts/ground_odts_app.h>
#include <lupnt/applications/lunar_sat_odts/satellite_odts_app.h>
#include <lupnt/applications/rover/surface_rover_nav_app.h>
#include <lupnt/simulations/simulation.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // <root>/data/LuPNT_data -> <root>/configs/<name>. Independent of the cwd.
  std::filesystem::path ConfigPath(const std::string& name) {
    return GetDataPath().parent_path().parent_path() / "configs" / name;
  }

  bool ConfigExists(const std::string& name) { return std::filesystem::exists(ConfigPath(name)); }

  // The south-pole DEM tile the rover/lander scenarios crop from. Present in the
  // repo data set, but guard so a fresh checkout without it skips cleanly.
  bool Site01DemPresent() {
    return std::filesystem::exists(GetDataPath() / "dem" / "LOLA_5mpp" / "Site01"
                                   / "Site01_final_adj_5mpp_surf.tif");
  }

  // True if every entry of a matrix/vector is finite (no NaN/Inf leaked through).
  bool AllFinite(const MatXd& m) {
    for (int i = 0; i < m.rows(); ++i)
      for (int j = 0; j < m.cols(); ++j)
        if (!std::isfinite(m(i, j))) return false;
    return true;
  }
}  // namespace

// ---------------------------------------------------------------------------
// ex9 / configs/ephemeris.yaml
// Exercises: simulations/world.cc, simulations/simulation.cc,
// applications/ephemeris/ephemeris_app.cc (+ ephemeris_basis / lunanet_*),
// agents/constellation + numerics integrators + environment force model.
// ---------------------------------------------------------------------------
TEST_CASE("simulations.integration.ephemeris", "[integration]") {
  if (!ConfigExists("ephemeris.yaml")) {
    SKIP("configs/ephemeris.yaml not found");
  }

  Config cfg = YAML::LoadFile(ConfigPath("ephemeris.yaml").string());

  // Shorten the truth arc + sweep so the fit study runs in seconds (in-memory only).
  Config app = cfg["agents"]["EphemerisManager"]["application"];
  app["duration_days"] = 0.1;  // ~2.4 h truth arc
  app["sample_dt_s"] = 60.0;
  app["fit_window_minutes"] = std::vector<double>{30.0};
  app["num_windows"] = 2;
  app["integration_step_s"] = 60.0;
  app["generate_almanac"] = true;

  Simulation sim(cfg);
  sim.Run();

  auto* eph = dynamic_cast<EphemerisApp*>(sim.GetAgent("EphemerisManager")->GetApplication().get());
  REQUIRE(eph != nullptr);

  const EphemerisResults& res = eph->GetResults();
  REQUIRE_FALSE(res.cartesian_results.empty());
  for (const EphemerisWindowResult& r : res.cartesian_results) {
    REQUIRE(r.num_params > 0);
    REQUIRE(r.total_bits > 0);
    REQUIRE(std::isfinite(r.pos_rms_m));
    REQUIRE(r.pos_rms_m >= 0.0);
    REQUIRE(std::isfinite(r.vel_rms_mps));
  }
}

// ---------------------------------------------------------------------------
// ex7 / configs/ground_station_odts.yaml
// Exercises: simulations/world.cc + simulation.cc, agents/satellite +
// ground_station, applications/ground_station/{tracking,manager}_app.cc,
// measurements (range/range-rate), filters/{batch_filter,srif}.cc.
// ---------------------------------------------------------------------------
TEST_CASE("simulations.integration.ground_station_odts", "[integration]") {
  if (!ConfigExists("ground_station_odts.yaml")) {
    SKIP("configs/ground_station_odts.yaml not found");
  }

  Config cfg = YAML::LoadFile(ConfigPath("ground_station_odts.yaml").string());
  cfg["duration"] = "03:00:00";  // short window of DSN passes (keeps the run to a few s)
  cfg["log_level"] = "WARNING";
  // Fewer batch iterations: this is a smoke/coverage run, not an accuracy study.
  cfg["agents"]["gs_manager"]["application"]["batch_max_iterations"] = 5;

  Simulation sim(cfg);
  sim.Run();

  auto mgr = std::dynamic_pointer_cast<GroundStationManagerApp>(
      sim.GetAgent("gs_manager")->GetApplication());
  REQUIRE(mgr != nullptr);
  REQUIRE(mgr->StationNames().size() == 3);

  // The DSN complexes should see the target at least once over 3 h, so the batch
  // filter aggregates measurements and produces an estimate + covariance.
  REQUIRE(mgr->NumMeasurements() > 0);
  REQUIRE(mgr->HasSolved());

  const Vec6d x0_est = mgr->X0Estimated();
  const Vec6d x0_true = mgr->X0True();
  REQUIRE(x0_est.allFinite());
  REQUIRE(x0_true.allFinite());

  const Mat6d cov = mgr->Covariance();
  REQUIRE(cov.allFinite());
  // A valid covariance has a non-negative position variance trace.
  REQUIRE(cov.topLeftCorner(3, 3).trace() > 0.0);

  // The estimate should not be wildly diverged from truth (loose sanity bound).
  const double pos_err = (x0_est.head(3) - x0_true.head(3)).norm();
  REQUIRE(std::isfinite(pos_err));
  REQUIRE(pos_err < 1.0e6);  // < 1000 km: the batch converged to a real orbit
}

// ---------------------------------------------------------------------------
// ex8 / configs/isl_odts_distributed.yaml
// Exercises: simulations/world.cc + simulation.cc, agents/spacecraft +
// surface_station, applications/lunar_sat_odts/{satellite_odts,ground_odts,
// station_beacon_sensor}.cc, devices/clock, measurements + Schmidt-EKF filters.
// ---------------------------------------------------------------------------
TEST_CASE("simulations.integration.isl_odts_distributed", "[integration]") {
  if (!ConfigExists("isl_odts_distributed.yaml")) {
    SKIP("configs/isl_odts_distributed.yaml not found");
  }

  Config cfg = YAML::LoadFile(ConfigPath("isl_odts_distributed.yaml").string());
  cfg["duration"] = "01:00:00";  // ~1 h -> a handful of filter epochs
  cfg["log_level"] = "WARNING";

  Simulation sim(cfg);
  sim.Run();

  const std::vector<std::string> sat_names = {"SV-1", "SV-2", "SV-3", "SV-4", "SV-5"};
  int N = -1;
  for (const std::string& sn : sat_names) {
    auto* app = dynamic_cast<SatelliteOdtsApp*>(sim.GetAgent(sn)->GetApplication().get());
    REQUIRE(app != nullptr);
    const MatXd& est = app->OwnEstimate();
    const MatXd& tru = app->TruthState();
    REQUIRE(est.rows() > 0);
    REQUIRE(est.rows() == tru.rows());
    REQUIRE(est.cols() == 8);  // [r(3) v(3) clock_bias clock_drift]
    REQUIRE(AllFinite(est));
    REQUIRE(AllFinite(tru));
    if (N < 0) N = static_cast<int>(est.rows());
    REQUIRE(static_cast<int>(est.rows()) == N);  // all onboard filters ran the same grid
  }
  REQUIRE(N > 0);

  // The centralized ground filter (station pseudoranges only) ran on the same grid.
  auto* ground = dynamic_cast<GroundOdtsApp*>(sim.GetAgent("gs_manager")->GetApplication().get());
  REQUIRE(ground != nullptr);
  const MatXd& estc = ground->EstCentral();
  REQUIRE(estc.rows() == N);
  REQUIRE(estc.cols() == 8 * static_cast<int>(sat_names.size()));
  REQUIRE(AllFinite(estc));
}

// ---------------------------------------------------------------------------
// ex10 / configs/surface_rover_nav.yaml
// Exercises: simulations/world.cc (DEM + gravity + ENU frame), agents/rover +
// constellation, applications/rover/surface_rover_nav_app.cc, devices/imu,
// measurements/surface_measurements, error-state INS EKF.
// Data-gated on the Site01 LOLA DEM tile.
// ---------------------------------------------------------------------------
TEST_CASE("simulations.integration.surface_rover_nav", "[integration]") {
  if (!ConfigExists("surface_rover_nav.yaml")) {
    SKIP("configs/surface_rover_nav.yaml not found");
  }
  if (!Site01DemPresent()) {
    SKIP("Site01 LOLA DEM tile not present (data-gated)");
  }

  Config cfg = YAML::LoadFile(ConfigPath("surface_rover_nav.yaml").string());
  cfg["duration"] = "00:01:00";  // 60 s
  cfg["log_level"] = "WARNING";
  Config app = cfg["agents"]["Rover"]["application"];
  app["duration_s"] = 60.0;  // match the shortened sim window
  app["dt_s"] = 1.0;

  Simulation sim(cfg);
  sim.Run();

  auto* rover = dynamic_cast<SurfaceRoverNavApp*>(sim.GetAgent("Rover")->GetApplication().get());
  REQUIRE(rover != nullptr);

  const VecXd& t = rover->time_series();
  REQUIRE(t.size() > 1);
  const VecXd& err = rover->pos_err_norm();
  REQUIRE(err.size() == t.size());
  REQUIRE(AllFinite(err));
  // The rover ran the full precomputed grid (duration_s / dt_s + 1 epochs).
  REQUIRE(t.size() == 61);
  REQUIRE(rover->pos_sigma_enu().rows() == t.size());
  REQUIRE(AllFinite(rover->pos_err_enu()));
}

// ---------------------------------------------------------------------------
// ex11 / configs/lander_nav.yaml
// Exercises: simulations/world.cc, agents/lander (multi-app), applications/
// lander/{lander_gnc_app,lander_nav_app}.cc, devices/imu, altimeter/crater/
// LunaNet measurements, error-state INS MEKF.
// Data-gated on the Site01 LOLA DEM tile.
// ---------------------------------------------------------------------------
TEST_CASE("simulations.integration.lander_nav", "[integration]") {
  if (!ConfigExists("lander_nav.yaml")) {
    SKIP("configs/lander_nav.yaml not found");
  }
  if (!Site01DemPresent()) {
    SKIP("Site01 LOLA DEM tile not present (data-gated)");
  }

  Config cfg = YAML::LoadFile(ConfigPath("lander_nav.yaml").string());
  cfg["duration"] = "00:00:30";  // 30 s of powered descent
  cfg["log_level"] = "WARNING";
  // The guidance app (index 0) owns the truth-trajectory duration the nav app reads.
  cfg["agents"]["Lander"]["applications"][0]["duration_s"] = 30.0;

  Simulation sim(cfg);
  sim.Run();

  auto nav = sim.GetAgent("Lander")->GetApplicationByName("LanderNavApp");
  auto* lander = dynamic_cast<LanderNavApp*>(nav.get());
  REQUIRE(lander != nullptr);

  const VecXd& t = lander->time_series();
  REQUIRE(t.size() > 1);
  const VecXd& err = lander->pos_err_norm();
  REQUIRE(err.size() == t.size());
  REQUIRE(AllFinite(err));
  REQUIRE(AllFinite(lander->pos_err_enu()));
  REQUIRE(lander->covariance().allFinite());
}
