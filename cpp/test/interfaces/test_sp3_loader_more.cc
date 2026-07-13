#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <filesystem>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// Deepens sp3_loader.cc coverage beyond test_sp3_loader_fixture.cc: several
// FilenameForEpoch/UrlForEpoch epochs (exercising the DayOfYear leap-year
// branch + the pre-GPS-epoch / pre-1962-week guard throws), the Chebyshev
// interpolation continuity / velocity-derivative path, multiple satellites,
// the boundary (t_min/t_max) queries, and the duplicate-epoch de-duplication
// merge path. Uses the committed trimmed SP3 fixture (no network).
namespace {
  std::filesystem::path FixtureFile() {
    return std::filesystem::path(LUPNT_TEST_FIXTURES_DIR) / "gnss" / "COD0MGXFIN_trimmed.SP3";
  }
}  // namespace

const double epsilon_more = 1e-6;

TEST_CASE("interfaces.sp3_fixture_more.download_helpers_epochs") {
  // Filenames are a pure function of the (GPS-time) calendar date, so they can
  // be cross-checked exactly. Midday epochs are used so the UTC->GPS offset
  // never crosses a day boundary. These cases exercise the leap-year branch of
  // the internal DayOfYear helper.

  // Non-leap year, last day of the year -> DOY 365.
  const Real t_2025_dec31 = GregorianToTime(2025, 12, 31, 12, 0, 0.0);
  REQUIRE(Sp3Loader::FilenameForEpoch(t_2025_dec31, Time::UTC)
          == "COD0MGXFIN_20253650000_01D_05M_ORB.SP3");

  // Leap year, Feb 29 -> DOY 60 (month == 2, leap-day itself, no +1 applied).
  const Real t_2024_feb29 = GregorianToTime(2024, 2, 29, 12, 0, 0.0);
  REQUIRE(Sp3Loader::FilenameForEpoch(t_2024_feb29, Time::UTC)
          == "COD0MGXFIN_20240600000_01D_05M_ORB.SP3");

  // Leap year, Mar 1 -> DOY 61 (leap && month > 2 adds one).
  const Real t_2024_mar01 = GregorianToTime(2024, 3, 1, 12, 0, 0.0);
  REQUIRE(Sp3Loader::FilenameForEpoch(t_2024_mar01, Time::UTC)
          == "COD0MGXFIN_20240610000_01D_05M_ORB.SP3");

  // Non-leap year, Mar 1 -> DOY 60.
  const Real t_2023_mar01 = GregorianToTime(2023, 3, 1, 12, 0, 0.0);
  REQUIRE(Sp3Loader::FilenameForEpoch(t_2023_mar01, Time::UTC)
          == "COD0MGXFIN_20230600000_01D_05M_ORB.SP3");

  // The URL is the modern CDDIS COD/MGEX layout: gnss/products/<week>/<file>.gz.
  const std::string url = Sp3Loader::UrlForEpoch(t_2025_dec31, Time::UTC);
  REQUIRE(url.rfind("https://cddis.nasa.gov/archive/gnss/products/", 0) == 0);
  REQUIRE(url.find("COD0MGXFIN_20253650000_01D_05M_ORB.SP3.gz") != std::string::npos);
  REQUIRE(url.substr(url.size() - 3) == ".gz");
}

TEST_CASE("interfaces.sp3_fixture_more.download_helpers_reject_old_epochs") {
  // Pre-GPS-epoch (before 1980-01-06) is rejected by the delta >= 0 guard.
  const Real t_pre_gps = GregorianToTime(1979, 1, 1, 12, 0, 0.0);
  REQUIRE_THROWS(Sp3Loader::FilenameForEpoch(t_pre_gps, Time::UTC));
  REQUIRE_THROWS(Sp3Loader::UrlForEpoch(t_pre_gps, Time::UTC));

  // Post-GPS but before GPS week 1962 (~Nov 2017) is rejected by the
  // "modern CDDIS COD/MGEX products only" guard.
  const Real t_old = GregorianToTime(2000, 1, 1, 12, 0, 0.0);
  REQUIRE_THROWS(Sp3Loader::FilenameForEpoch(t_old, Time::UTC));
  REQUIRE_THROWS(Sp3Loader::UrlForEpoch(t_old, Time::UTC));
}

TEST_CASE("interfaces.sp3_fixture_more.interpolation_continuity_and_derivative") {
  Sp3Loader loader(FixtureFile());
  const std::string sat_id = "G01";
  auto [t_min, t_max] = loader.GetTimeSpan(sat_id);

  // Boundary queries are inclusive (t_min and t_max both valid).
  Vec6 rv_min = loader.GetPosVel(sat_id, Real(t_min));
  Vec6 rv_max = loader.GetPosVel(sat_id, Real(t_max));
  REQUIRE(std::isfinite(rv_min.head(3).norm().val()));
  REQUIRE(std::isfinite(rv_max.head(3).norm().val()));

  // The Chebyshev-fit velocity must match a central finite difference of the
  // interpolated position (validates the analytic-derivative velocity path).
  Real t_mid = Real(0.5 * (t_min + t_max));
  const double dt = 1.0;
  Vec6 rv = loader.GetPosVel(sat_id, t_mid);
  Vec6 rv_p = loader.GetPosVel(sat_id, Real(t_mid.val() + dt));
  Vec6 rv_m = loader.GetPosVel(sat_id, Real(t_mid.val() - dt));
  Vec3 v_fd = (rv_p.head(3) - rv_m.head(3)) / (2.0 * dt);
  REQUIRE((v_fd - rv.tail(3)).norm().val() < 1.0);  // < 1 m/s

  // Position must be spatially continuous over the sample step.
  REQUIRE((rv_p.head(3) - rv.head(3)).norm().val() < 1.0e4);  // < 10 km over 1 s
}

TEST_CASE("interfaces.sp3_fixture_more.multiple_satellites") {
  Sp3Loader loader(FixtureFile());

  // GPS (G02/G03) and Galileo (E02) all parse and interpolate to plausible MEO
  // radii / speeds. This exercises the per-satellite table / model machinery
  // for satellites other than the single G01 covered elsewhere.
  for (const std::string& sat : {std::string("G02"), std::string("G03"), std::string("E02")}) {
    REQUIRE(loader.HasSatellite(sat));
    auto [t_min, t_max] = loader.GetTimeSpan(sat);
    REQUIRE(t_max > t_min);

    Vec6 rv = loader.GetPosVel(sat, Real(0.5 * (t_min + t_max)));
    Real r = rv.head(3).norm();
    Real v = rv.tail(3).norm();
    REQUIRE(r.val() > 1.8e7);
    REQUIRE(r.val() < 3.5e7);
    REQUIRE(v.val() > 1.0e3);
    REQUIRE(v.val() < 1.0e4);
  }
}

TEST_CASE("interfaces.sp3_fixture_more.duplicate_epoch_dedup_merge") {
  const std::filesystem::path sp3_file = FixtureFile();
  const std::string sat_id = "G01";

  // Loading the *same* file twice makes every epoch a duplicate; the de-dup
  // logic in LoadFile must collapse them back so the span and the interpolated
  // state are identical to the single-file load (no doubling / corruption).
  Sp3Loader loader_single(sp3_file);
  auto [s_min, s_max] = loader_single.GetTimeSpan(sat_id);

  Sp3Loader loader_dup(std::vector<std::filesystem::path>{sp3_file, sp3_file});
  auto [d_min, d_max] = loader_dup.GetTimeSpan(sat_id);
  REQUIRE_THAT(d_min, WithinAbs(s_min, 1e-6));
  REQUIRE_THAT(d_max, WithinAbs(s_max, 1e-6));

  Real t_mid = Real(0.5 * (s_min + s_max));
  Vec6 rv_single = loader_single.GetPosVel(sat_id, t_mid);
  Vec6 rv_dup = loader_dup.GetPosVel(sat_id, t_mid);
  RequireNear(rv_dup, rv_single, epsilon_more);

  // Incremental LoadFile onto an already-loaded file matches the vector ctor.
  Sp3Loader loader_incr(sp3_file);
  loader_incr.LoadFile(sp3_file);
  auto [i_min, i_max] = loader_incr.GetTimeSpan(sat_id);
  REQUIRE_THAT(i_min, WithinAbs(s_min, 1e-6));
  REQUIRE_THAT(i_max, WithinAbs(s_max, 1e-6));
}
