#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "utils.cc"
using namespace lupnt;
using namespace Catch::Matchers;

const Real M_ref = 0.06981317008;
const Real e_ref = 0.72000000000;
const Real E_ref = 0.24318719637;
const Real nu_ref = Mean2TrueAnomaly(M_ref, e_ref);

TEST_CASE("OrbitState", "Anomaly") {
  const double eps = 1e-9;

  // Example 2-2 (Kepler’s equation)
  const Real M_ref = 0.06981317008;
  const Real e_ref = 0.72000000000;
  const Real E_ref = 0.24318719637;
  const Real nu_ref = Mean2TrueAnomaly(M_ref, e_ref);

  RequireNearReal(Mean2EccAnomaly(M_ref, e_ref), E_ref, eps);
  RequireNearReal(Ecc2MeanAnomaly(E_ref, e_ref), M_ref, eps);
  RequireNearReal(Ecc2TrueAnomaly(E_ref, e_ref), nu_ref, eps);
  RequireNearReal(True2EccAnomaly(nu_ref, e_ref), E_ref, eps);
  RequireNearReal(Mean2TrueAnomaly(M_ref, e_ref), nu_ref, eps);
  RequireNearReal(True2MeanAnomaly(nu_ref, e_ref), M_ref, eps);
}

TEST_CASE("OrbitState", "Conversions") {
  const double eps = 1e-3;
  const double eps_ecc = 1e-3;

  // Example 2-3 (Osculating Elements)
  Vec6 rv_ref(10e3, 40e3, -5e3, -1.5, 1, -0.1);  // [km, km/s]
  Vec6 coe_ref(25015.181, 0.7079772, 6.971, 173.290, 91.553, 144.225);
  coe_ref.segment(2, 4) *= RAD;

  Vec6 rv = Classical2Cart(coe_ref, GM_EARTH);
  Vec6 coe = Cart2Classical(rv, GM_EARTH);

  RequireNearReal(coe(1), coe_ref(1), eps_ecc);
  RequireNearRealVec(coe, coe_ref, eps);
  RequireNearRealVec(rv, rv_ref, eps);
}
