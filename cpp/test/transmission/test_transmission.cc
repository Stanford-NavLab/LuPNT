#include <lupnt/transmission/transmission.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("transmission.transmission") {
  auto data = MakePtr<SimpleGnssData>();
  data->rv_tx = Cart6(Vec6(1.0, 2.0, 3.0, 0.1, 0.2, 0.3), Frame::GCRF);
  data->clk_tx = ClockState2(Vec2(1.0e-6, 2.0e-9));

  SimpleGnssTransmission transmission;
  transmission.data = data;
  transmission.P_tx = 42.0;

  REQUIRE(transmission.data == data);
  REQUIRE_THAT(transmission.P_tx.val(), Catch::Matchers::WithinAbs(42.0, epsilon));
  REQUIRE_THAT(transmission.data->rv_tx(0).val(), Catch::Matchers::WithinAbs(1.0, epsilon));
  REQUIRE_THAT(transmission.data->clk_tx.b().val(), Catch::Matchers::WithinAbs(1.0e-6, 1e-12));
}
