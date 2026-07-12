// Unit tests for `ParseForceModelSpec` -- the shared parser that turns a unified `force_model:`
// YAML block (bodies list + relativity + SRP) into scalar `ForceModelSpec` fields, used so every
// scenario specifies dynamics the same way.
#include <lupnt/dynamics/numerical_orbit_dynamics.h>

#include <catch2/catch_test_macros.hpp>
#include <yaml-cpp/yaml.h>

#include "../utils.cc"

using namespace lupnt;

TEST_CASE("dynamics.force_model_spec") {
  SECTION("full bodies list + relativity + SRP") {
    YAML::Node n = YAML::Load(R"(
bodies:
  - MOON: { n: 20, m: 18 }
  - EARTH: {}
  - SUN: {}
relativity: true
CR: 1.3
area: 0.02
mass: 5.0
)");
    ForceModelSpec s = ParseForceModelSpec(n);
    REQUIRE(s.moon_degree == 20);
    REQUIRE(s.moon_order == 18);
    REQUIRE(s.include_earth);
    REQUIRE(s.include_sun);
    REQUIRE(s.relativity);
    REQUIRE(s.has_srp);
    REQUIRE(s.srp_cr == 1.3);
    REQUIRE(s.srp_area_m2 == 0.02);
    REQUIRE(s.srp_mass_kg == 5.0);
  }

  SECTION("point-mass Moon only -> degree 0, no third bodies, no SRP") {
    YAML::Node n = YAML::Load("bodies: [ { MOON: {} } ]");
    ForceModelSpec s = ParseForceModelSpec(n);
    REQUIRE(s.moon_degree == 0);
    REQUIRE(s.moon_order == 0);
    REQUIRE_FALSE(s.include_earth);
    REQUIRE_FALSE(s.include_sun);
    REQUIRE_FALSE(s.relativity);
    REQUIRE_FALSE(s.has_srp);
  }

  SECTION("`use_relativity` spelling is also accepted") {
    YAML::Node n = YAML::Load("bodies: [ {EARTH: {}} ]\nuse_relativity: true");
    ForceModelSpec s = ParseForceModelSpec(n);
    REQUIRE(s.relativity);
    REQUIRE(s.include_earth);
    REQUIRE(s.moon_degree == 0);  // no MOON entry
  }

  SECTION("empty / missing block -> all defaults off") {
    YAML::Node n = YAML::Load("{}");
    ForceModelSpec s = ParseForceModelSpec(n);
    REQUIRE(s.moon_degree == 0);
    REQUIRE_FALSE(s.include_earth);
    REQUIRE_FALSE(s.include_sun);
    REQUIRE_FALSE(s.relativity);
    REQUIRE_FALSE(s.has_srp);
  }

  SECTION("CR present without area/mass -> SRP flagged, area/mass keep defaults") {
    YAML::Node n = YAML::Load("bodies: [ {MOON: {n: 2, m: 0}} ]\nCR: 1.0");
    ForceModelSpec s = ParseForceModelSpec(n);
    REQUIRE(s.moon_degree == 2);
    REQUIRE(s.moon_order == 0);
    REQUIRE(s.has_srp);
    REQUIRE(s.srp_cr == 1.0);
  }
}
