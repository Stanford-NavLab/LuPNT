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

TEST_CASE("measurements.comms_utils.visibility_geometry") {
  const Real R = Real(R_EARTH);

  // Two elevated points on the same side of the body: clear line of sight.
  Vec3 a(30000e3, 0.0, 0.0);
  Vec3 b(30000e3, 1000e3, 0.0);
  REQUIRE(ComputeVisibility(a, b, R));

  // Line of sight through the body's center is occluded.
  Vec3 c(-30000e3, 0.0, 0.0);
  REQUIRE_FALSE(ComputeVisibility(a, c, R));

  // Occluding-body offset from the origin: same geometry, translated.
  Vec3 r_body(100e3, 200e3, -50e3);
  REQUIRE(ComputeVisibility(Vec3(a + r_body), Vec3(b + r_body), R, r_body));
  REQUIRE_FALSE(ComputeVisibility(Vec3(a + r_body), Vec3(c + r_body), R, r_body));
}

TEST_CASE("measurements.comms_utils.visibility_elevation_mask") {
  const Real R = Real(R_EARTH);
  Vec3 surface(R_EARTH, 0.0, 0.0);

  // Satellite along the local horizontal (~0 deg elevation) falls below the mask.
  Vec3 horizon_sat(R_EARTH, 30000e3, 0.0);
  REQUIRE_FALSE(ComputeVisibility(surface, horizon_sat, R));
  // Symmetric handling when the surface point is the second argument.
  REQUIRE_FALSE(ComputeVisibility(horizon_sat, surface, R));

  // Satellite well above the horizon (~45 deg) is visible.
  Vec3 high_sat(R_EARTH + 20000e3, 20000e3, 0.0);
  REQUIRE(ComputeVisibility(surface, high_sat, R));

  // The elevation mask is configurable: the same ~45 deg satellite passes a 1 deg
  // mask but is rejected by an 80 deg mask.
  REQUIRE(ComputeVisibility(surface, high_sat, R, Vec3::Zero(), 10e3, 1.0));
  REQUIRE_FALSE(ComputeVisibility(surface, high_sat, R, Vec3::Zero(), 10e3, 80.0));
}
