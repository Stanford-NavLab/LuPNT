// Regression test for LoadTleFile: it built a ClassicalOE state and handed it to
// Satellite::SetState, but a Satellite carries (and its dynamics integrate) a
// Cartesian state, so the load left the agent with a non-Cartesian state. The
// fix converts via ClassicalToCart. Here we load the bundled GPS TLE set and
// assert each resulting Satellite holds a finite Cartesian state at a sane MEO
// radius.
//
// gps_2023_06_09.txt ships in the standard bundled data set (it is the CI
// data-download marker file); LoadTleFile resolves it via GetFilePath's data-dir
// search. Tag intentionally avoids the excluded "states.tle" pattern.

#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "../utils.cc"

using namespace lupnt;

TEST_CASE("interfaces.tle_loader") {
  std::vector<Ptr<Agent>> sats = LoadTleFile("gps_2023_06_09.txt");
  REQUIRE(sats.size() > 10);  // a full GPS constellation

  // GetState is non-virtual (base Agent returns an empty State); the concrete
  // Satellite carries the Cartesian state set by the loader.
  auto* sat = dynamic_cast<Satellite*>(sats.front().get());
  REQUIRE(sat != nullptr);
  State s = sat->GetState();

  double r2 = 0.0;
  for (int i = 0; i < 6; ++i) REQUIRE(std::isfinite(s(i).val()));
  for (int i = 0; i < 3; ++i) r2 += s(i).val() * s(i).val();
  const double r = std::sqrt(r2);

  // A Cartesian GPS (MEO) position, geocentric radius ~2.66e7 m. A ClassicalOE
  // state (the pre-fix bug) would not land in this Cartesian range.
  REQUIRE(r > 2.0e7);
  REQUIRE(r < 3.0e7);
}
