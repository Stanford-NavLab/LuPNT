#include <lupnt/interfaces/yaml.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("interfaces.yaml") {
  Config config = YAML::Load(R"(
real: 1.25
vec: [1.0, 2.0, 3.0]
mat:
  - [1.0, 2.0]
  - [3.0, 4.0]
)");
  Real real = config["real"].as<Real>();
  Vec3d vec = config["vec"].as<Vec3d>();
  Mat2d mat = config["mat"].as<Mat2d>();

  REQUIRE_THAT(real.val(), Catch::Matchers::WithinAbs(1.25, epsilon));
  REQUIRE(vec.isApprox(Vec3d(1.0, 2.0, 3.0), epsilon));
  REQUIRE(mat.isApprox((Mat2d() << 1.0, 2.0, 3.0, 4.0).finished(), epsilon));

  YAML::Node encoded = YAML::convert<Vec3d>::encode(vec);
  REQUIRE(encoded.IsSequence());
  REQUIRE(encoded.size() == 3);
}
