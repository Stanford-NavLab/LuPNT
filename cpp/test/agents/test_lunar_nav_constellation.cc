#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // Shared point-mass Moon dynamics block every satellite inherits.
  const char* kDynamics
      = "dynamics:\n"
        "  class: NBodyDynamics\n"
        "  frame: MOON_CI\n"
        "  integrator: RKF45\n"
        "  dt: 60.0\n"
        "  bodies:\n"
        "    - MOON: {}\n";
}  // namespace

TEST_CASE("agents.lunar_nav_constellation.explicit_list") {
  std::string yaml = std::string("name: LuNaNet\n") + kDynamics
                     + "satellites:\n"
                       "  - { name: SV-1, r0_m: [6500000.0, 0.0, 0.0], v0_mps: [0.0, 1600.0, 0.0] }\n"
                       "  - { name: SV-2, r0_m: [0.0, 6500000.0, 0.0], v0_mps: [-1600.0, 0.0, 0.0] }\n"
                       "  - { name: SV-3, r0_m: [0.0, 0.0, 6500000.0], v0_mps: [1600.0, 0.0, 0.0] }\n";
  Config config = YAML::Load(yaml);
  LunarNavConstellation constellation(config);

  REQUIRE(constellation.NumSatellites() == 3);
  REQUIRE(constellation.SatelliteName(0) == "SV-1");
  REQUIRE(constellation.SatelliteName(1) == "SV-2");
  REQUIRE(constellation.SatelliteName(2) == "SV-3");

  // At t = 0 the truth state equals the configured initial state (no propagation).
  Cart6 s0 = constellation.GetSatelliteStateAt(0, 0.0);
  REQUIRE_THAT(s0.r()(0).val(), WithinAbs(6500000.0, 1.0e-6));
  REQUIRE_THAT(s0.r()(1).val(), WithinAbs(0.0, 1.0e-6));
  REQUIRE_THAT(s0.v()(1).val(), WithinAbs(1600.0, 1.0e-6));

  Cart6 s1 = constellation.GetSatelliteStateAt(1, 0.0);
  REQUIRE_THAT(s1.r()(1).val(), WithinAbs(6500000.0, 1.0e-6));
  REQUIRE_THAT(s1.v()(0).val(), WithinAbs(-1600.0, 1.0e-6));

  // GetStateAt (constellation-level) returns the first satellite's state.
  Cart6 s_default = constellation.GetStateAt(0.0);
  REQUIRE_THAT((s_default.r() - s0.r()).norm().val(), WithinAbs(0.0, 1.0e-6));
}

TEST_CASE("agents.lunar_nav_constellation.walker") {
  // 2 planes x 3 satellites = 6, circular (e=0) frozen orbit; |r| == a for every satellite.
  const double a = 6541000.0;
  std::string yaml = std::string("name: Walker\n") + kDynamics
                     + "walker:\n"
                       "  n_planes: 2\n"
                       "  sats_per_plane: 3\n"
                       "  a: 6541000.0\n"
                       "  e: 0.0\n"
                       "  i: 56.0\n"
                       "  omega: 0.0\n"
                       "  raan0_deg: 0.0\n"
                       "  m0_deg: 0.0\n"
                       "  frame: MOON_CI\n"
                       "  name_prefix: NAV\n";
  Config config = YAML::Load(yaml);
  LunarNavConstellation constellation(config);

  REQUIRE(constellation.NumSatellites() == 6);
  REQUIRE(constellation.SatelliteName(0) == "NAV-1");
  REQUIRE(constellation.SatelliteName(5) == "NAV-6");

  // Each satellite sits on the circular orbit of radius a.
  for (int j = 0; j < constellation.NumSatellites(); ++j) {
    Cart6 s = constellation.GetSatelliteStateAt(j, 0.0);
    REQUIRE_THAT(s.r().norm().val(), WithinAbs(a, 1.0e-3));
  }

  // The three satellites in the first plane are spread 120 deg in mean/true anomaly, so their
  // position unit vectors differ.
  Cart6 sa = constellation.GetSatelliteStateAt(0, 0.0);
  Cart6 sb = constellation.GetSatelliteStateAt(1, 0.0);
  REQUIRE((sa.r() - sb.r()).norm().val() > 1.0e5);
}

TEST_CASE("agents.lunar_nav_constellation.invalid_config_throws") {
  // Missing both `satellites:` and `walker:` blocks.
  std::string yaml = std::string("name: Bad\n") + kDynamics;
  Config config = YAML::Load(yaml);
  REQUIRE_THROWS(LunarNavConstellation(config));
}
