// Engine-level coverage for the lunar-GNSS ODTS Monte-Carlo pipeline
// (`simulations/lunar_gnss_odts/lunar_gnss_odts_simulation.cc`). That engine is
// normally data-gated on downloaded SP3/BRDC/ANTEX products plus a plasma model,
// so it gets ~0 CI coverage. This test drives the full public struct API --
// ResolveLunarGnssODTSConfigForRun -> PrecomputeLunarGnssODTSLinks (Stage 1 link
// geometry) -> RunLunarGnssODTSMonteCarlo (the UDU-EKF) -- using ONLY committed
// fixtures and no network:
//   * a YUMA almanac (`current.alm`) seeds the GPS constellation, which is then
//     numerically propagated (no date-specific SP3/BRDC needed), and
//   * a trimmed ANTEX (`igs20_trimmed.atx`, GPS-only) supplies antenna PCOs.
// Plasma truth/filter are disabled (no plasma data). A short 300 s arc keeps the
// run to a couple of seconds. This is a code-path/coverage smoke test, not an
// accuracy study, so we only assert the summary fields are finite and
// non-negative rather than any convergence bound.
//
// `LUPNT_TEST_FIXTURES_DIR` is defined by `cpp/test/CMakeLists.txt`.

#include <lupnt/lupnt.h>
#include <lupnt/simulations/lunar_gnss_odts/lunar_gnss_odts_simulation.h>

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <filesystem>

#include "../utils.cc"

using namespace lupnt;

namespace {
  std::filesystem::path GnssFixture(const std::string& name) {
    return std::filesystem::path(LUPNT_TEST_FIXTURES_DIR) / "gnss" / name;
  }
}  // namespace

TEST_CASE("simulations.lunar_gnss_odts") {
  const std::filesystem::path almanac_file = GnssFixture("current.alm");
  const std::filesystem::path antex_file = GnssFixture("igs20_trimmed.atx");
  REQUIRE(std::filesystem::exists(almanac_file));
  REQUIRE(std::filesystem::exists(antex_file));

  // Fixed, unique scratch dir. remove_all up front so reruns are deterministic
  // (no stale/mismatched link cache is left behind between invocations).
  const std::filesystem::path out_dir
      = std::filesystem::temp_directory_path() / "lupnt_gnss_odts_engine_test";
  std::error_code ec;
  std::filesystem::remove_all(out_dir, ec);

  LunarGnssODTSConfig cfg;
  cfg.seed = 42;
  cfg.monte_carlo_runs = 1;

  // Short 300 s arc (receiver app runs at its default 1 Hz), so the whole
  // pipeline -- precompute + EKF over the epochs -- runs in a couple of seconds.
  cfg.duration_s = 300.0;
  cfg.dt_s = 60.0;
  cfg.ephemeris_dt_s = 60.0;
  cfg.start_epoch_utc = "2026-01-14T00:00:00";  // matches the ANTEX fixture validity

  cfg.output_dir = out_dir;
  cfg.links_file = out_dir / "precomputed_links.csv";
  cfg.delays_file = out_dir / "precomputed_delays.csv";

  // Almanac-seeded, numerically-propagated GPS constellation: no date-specific
  // SP3/BRDC required. The ANTEX fixture is GPS-only, so keep Galileo off.
  cfg.constellation.source = "almanac";
  cfg.constellation.almanac_file = almanac_file;
  cfg.constellation.antex_file = antex_file;
  cfg.constellation.use_all_gps = true;
  cfg.constellation.include_galileo = false;

  // No plasma data available -> disable ionosphere/plasmasphere truth and filter model.
  cfg.plasma.simulate_truth = false;
  cfg.plasma.model_in_filter = false;

  // No SRP; modest (but stable) gravity for speed. Keep the RK4 step at its
  // default (large steps destabilize MEO propagation).
  cfg.use_srp_truth = false;
  cfg.use_srp_filter = false;
  cfg.moon_gravity_degree_truth = 4;
  cfg.moon_gravity_order_truth = 4;
  cfg.moon_gravity_degree_filter = 4;
  cfg.moon_gravity_order_filter = 4;
  cfg.moon_gravity_degree_constellation = 4;
  cfg.moon_gravity_order_constellation = 4;

  // Guarantee measurements are formed regardless of C/N0 geometry.
  cfg.design.apply_cn0_threshold = false;

  // Resolve (no-op for the almanac source) then precompute the link cache that
  // RunLunarGnssODTSMonteCarlo requires (it throws if the cache is missing).
  ResolveLunarGnssODTSConfigForRun(cfg);
  PrecomputeLunarGnssODTSLinks(cfg);
  REQUIRE(std::filesystem::exists(cfg.links_file));

  std::vector<LunarGnssODTSSummary> summaries = RunLunarGnssODTSMonteCarlo(cfg);

  REQUIRE_FALSE(summaries.empty());
  const LunarGnssODTSSummary& s = summaries[0];
  REQUIRE(s.num_epochs > 0);

  // Coverage assertions only: the EKF error metrics are finite and non-negative.
  // This exercises the code path; do NOT assert accuracy/convergence.
  REQUIRE(std::isfinite(s.final_position_error_m));
  REQUIRE(s.final_position_error_m >= 0.0);
  REQUIRE(std::isfinite(s.rms_position_error_m));
  REQUIRE(s.rms_position_error_m >= 0.0);

  std::filesystem::remove_all(out_dir, ec);
}
