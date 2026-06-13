#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // The downloaded SP3/BRDC GNSS files live under `<project_root>/output/gnss_files`
  // (as produced by `pylupnt.interfaces.gnss_file_loader.SP3Loader`/`BRDCLoader`,
  // see `python/pylupnt/interfaces/gnss_file_loader.py`), which does **not**
  // coincide with `GetOutputDir("gnss_files")` (== `GetDataPath()/output/gnss_files`)
  // when `$LUPNT_OUTPUT_PATH` is unset. `GetDataPath()` resolves to
  // `<project_root>/data/LuPNT_data`, so two `parent_path()` hops recover
  // `<project_root>`.
  std::filesystem::path GnssFilesDir() {
    return GetDataPath().parent_path().parent_path() / "output" / "gnss_files";
  }
}  // namespace

const double epsilon = 1e-6;

TEST_CASE("interfaces.sp3_loader") {
  const std::filesystem::path sp3_dir = GnssFilesDir() / "sp3";
  const std::filesystem::path sp3_file_1 = sp3_dir / "COD0MGXFIN_20260140000_01D_05M_ORB.SP3";
  const std::filesystem::path sp3_file_2 = sp3_dir / "COD0MGXFIN_20260150000_01D_05M_ORB.SP3";

  const std::string sat_id = "G01";  // GPS PRN 01 (present in both files, see SP3 header)

  SECTION("Default-constructed loader has no satellites loaded") {
    Sp3Loader loader;
    REQUIRE(loader.GetSatellites().empty());
    REQUIRE_FALSE(loader.HasSatellite(sat_id));
  }

  SECTION("Download helpers map UTC epochs to modern CDDIS COD/MGEX products") {
    const Real t_utc = GregorianToTime(2025, 1, 1, 0, 0, 0.0);
    REQUIRE(Sp3Loader::FilenameForEpoch(t_utc, Time::UTC)
            == "COD0MGXFIN_20250010000_01D_05M_ORB.SP3");
    REQUIRE(Sp3Loader::UrlForEpoch(t_utc, Time::UTC)
            == "https://cddis.nasa.gov/archive/gnss/products/2347/"
               "COD0MGXFIN_20250010000_01D_05M_ORB.SP3.gz");

    const Real t_midday_utc = GregorianToTime(2025, 1, 3, 12, 0, 0.0);
    REQUIRE(Sp3Loader::FilenameForEpoch(t_midday_utc, Time::UTC)
            == "COD0MGXFIN_20250030000_01D_05M_ORB.SP3");
  }

  SECTION("Single-file constructor parses satellites & ephemeris samples") {
    Sp3Loader loader(sp3_file_1);

    REQUIRE_FALSE(loader.GetSatellites().empty());
    REQUIRE(loader.HasSatellite(sat_id));
    REQUIRE_FALSE(loader.HasSatellite("X99"));

    auto [t_min, t_max] = loader.GetTimeSpan(sat_id);
    REQUIRE(t_max > t_min);

    Real t_mid = Real(0.5 * (t_min + t_max));

    Vec6 rv_ecef;
    Real clock_bias_s;
    loader.GetPosVelClock(sat_id, t_mid, rv_ecef, clock_bias_s);

    // GPS satellites orbit at an altitude of ~20,200 km -> geocentric radius ~26,560 km
    Real r = rv_ecef.head(3).norm();
    REQUIRE(r.val() > 2.0e7);
    REQUIRE(r.val() < 3.0e7);

    // GPS orbital speed is ~3.9 km/s
    Real v = rv_ecef.tail(3).norm();
    REQUIRE(v.val() > 1.0e3);
    REQUIRE(v.val() < 1.0e4);

    // SP3 clock products are accurate to well under a millisecond
    REQUIRE(std::abs(clock_bias_s.val()) < 1e-3);

    // GetPosVel returns the same position/velocity as GetPosVelClock
    Vec6 rv_only = loader.GetPosVel(sat_id, t_mid);
    RequireNear(rv_only, rv_ecef, epsilon);
  }

  SECTION("LoadFile / multi-file constructor merge & sort samples across days") {
    Sp3Loader loader_a(sp3_file_1);
    auto [t1_min, t1_max] = loader_a.GetTimeSpan(sat_id);

    loader_a.LoadFile(sp3_file_2);
    auto [t2_min, t2_max] = loader_a.GetTimeSpan(sat_id);

    // Loading the (consecutive) next day extends the covered time span
    REQUIRE(t2_min <= t1_min);
    REQUIRE(t2_max > t1_max);

    // The multi-file constructor produces the same merged span directly
    Sp3Loader loader_b(std::vector<std::filesystem::path>{sp3_file_1, sp3_file_2});
    auto [t3_min, t3_max] = loader_b.GetTimeSpan(sat_id);
    REQUIRE_THAT(t3_min, WithinAbs(t2_min, 1e-6));
    REQUIRE_THAT(t3_max, WithinAbs(t2_max, 1e-6));
  }

  SECTION("Queries for unknown satellites or out-of-span epochs throw") {
    Sp3Loader loader(sp3_file_1);
    Vec6 rv_ecef;
    Real clock_bias_s;

    REQUIRE_THROWS(loader.GetTimeSpan("X99"));
    REQUIRE_THROWS(loader.GetPosVelClock("X99", Real(0.0), rv_ecef, clock_bias_s));

    auto [t_min, t_max] = loader.GetTimeSpan(sat_id);
    REQUIRE_THROWS(loader.GetPosVelClock(sat_id, Real(t_min - 1.0e6), rv_ecef, clock_bias_s));
    REQUIRE_THROWS(loader.GetPosVelClock(sat_id, Real(t_max + 1.0e6), rv_ecef, clock_bias_s));
  }
}
