#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../data.cc"
#include "../utils.cc"
using namespace lupnt;
using namespace Catch::Matchers;

TEST_CASE("Lunar Mean Osculating") {
  Vec6 coe_mean = GetClassicalOE();
  Real J2 = J2_MOON;
  Real GM = GM_MOON;
  Vec6 coe_osc = Mean2Osculating(coe_mean, GM, J2);
  Vec6 coe_mean2 = Osculating2Mean(coe_osc, GM, J2);

  ClassicalOE coe_mean_state(coe_mean);
  ClassicalOE coe_osc_state = Mean2Osculating(coe_mean_state, GM, J2);
  ClassicalOE coe_mean2_state = Osculating2Mean(coe_osc_state, GM, J2);
}
