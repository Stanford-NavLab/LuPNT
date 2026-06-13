#include <lupnt/interfaces/cesium.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("interfaces.cesium") {
  VecX times(2);
  times << 0.0, 10.0;
  MatX3 positions(2, 3);
  positions << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

  Entity entity{"id",      "name",        times,        positions,
                {1, 2, 3}, "description", BodyId::MOON, "2024-01-01T00:00:00Z",
                7};
  REQUIRE(entity.id == "id");
  REQUIRE(entity.name == "name");
  REQUIRE(entity.times.size() == 2);
  REQUIRE(entity.positions.rows() == 2);
  REQUIRE(entity.color[2] == 3);
  REQUIRE(entity.body_id == BodyId::MOON);
  REQUIRE(entity.size == 7);
}
