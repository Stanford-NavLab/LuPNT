#include <lupnt/measurements/comms_utils.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("measurements.comms_utils") {
  DllParams dll{1.0, 0.02, 0.5, 1.0e-3, 2.0e6};
  PllParams pll{15.0, 0.02};
  FllParams fll{10.0, 0.02, 25.0};
  Real cn0 = 1.0e4;

  REQUIRE_THAT(SigmaDll(dll, cn0).val(), Catch::Matchers::WithinAbs(0.00501663898, 1e-10));
  REQUIRE_THAT(SigmaPll(pll, cn0).val(), Catch::Matchers::WithinAbs(0.03877821553, 1e-10));
  REQUIRE(SigmaFll(fll, cn0) > 0.0);
  REQUIRE(SigmaFll(fll, Real(1.0)) > SigmaFll(fll, cn0));

  ArrX cn0s = ArrX::Constant(2, 1, cn0);
  REQUIRE_THAT(SigmaDll(dll, cn0s)(0, 0).val(),
               Catch::Matchers::WithinAbs(SigmaDll(dll, cn0).val(), epsilon));
  REQUIRE_THAT(SigmaPll(pll, cn0s)(1, 0).val(),
               Catch::Matchers::WithinAbs(SigmaPll(pll, cn0).val(), epsilon));

  REQUIRE_THAT(FreeSpacePathLoss(Real(2.0), Real(C / (4.0 * PI))).val(),
               Catch::Matchers::WithinAbs(20.0 * std::log10(2.0), epsilon));
  REQUIRE_THAT(ParabolicAntennaGain(Real(2.0), Real(4.0), Real(12.0)).val(),
               Catch::Matchers::WithinAbs(9.0, epsilon));
}
