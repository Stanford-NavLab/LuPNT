#include <lupnt/conversions/time_conversions.h>
#include <lupnt/data/kernels.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// ============================================================================
// Cross-validation against Orekit
//
// The reference values below were generated once via Orekit and checked
// into `cpp/test/orekit/data/orekit_reference.json` by
// `cpp/test/orekit/gen_orekit_reference.py`. This test does **not**
// require Orekit/Java to run -- see `cpp/test/orekit/README.md` for how
// the fixture was produced and how to regenerate it.
//
// Both LuPNT and Orekit evaluate the JPL DE440 ephemerides, but through
// completely independent code paths (LuPNT's Chebyshev-coefficient reader
// vs Orekit's JPL-DE loader), so this test catches indexing, units,
// time-argument, and interpolation errors in LuPNT's ephemeris evaluation.
// ============================================================================

TEST_CASE("data.ephemeris_orekit_reference") {
  nlohmann::json data = LoadTestJson("orekit/data/orekit_reference.json");

  for (const auto& tc : data["ephemerides"]) {
    std::string epoch = tc["epoch_utc"].get<std::string>();
    std::string target_name = tc["target"].get<std::string>();
    Real t_tdb = ConvertTime(GregorianToTime(epoch), Time::UTC, Time::TDB);

    BodyId target = (target_name == "MOON") ? BodyId::MOON : BodyId::SUN;

    DYNAMIC_SECTION("epoch " << epoch << ", target " << target_name) {
      Vec6 rv = GetBodyPosVel(t_tdb, BodyId::EARTH, target, Frame::GCRF);

      // Both sides evaluate DE440-series data, but through separately
      // distributed copies (LuPNT's ASCII Chebyshev coefficients vs the
      // orekit-data bundle), which agree to ~5e-11 relative: observed
      // differences are ~0.1-0.35 m for the Moon (|r| ~ 4e8 m) and
      // ~1-7 m for the Sun (|r| ~ 1.5e11 m). The position tolerance
      // therefore scales with |r| (1e-10 * |r|, with a 1 m floor), which
      // still catches any real error -- wrong segment, km vs m, TDB vs TT
      // time argument -- as those show up at km scale or more.
      double r_norm = rv.head<3>().norm().val();
      double r_tol = std::max(1.0, 1e-10 * r_norm);
      for (int i = 0; i < 3; i++) {
        RequireNear(rv(i), Real(tc["r_gcrf"].at(i).get<double>()), r_tol);     // [m]
        RequireNear(rv(i + 3), Real(tc["v_gcrf"].at(i).get<double>()), 1e-5);  // [m/s]
      }
    }
  }
}
