#pragma once

#include "orbit_states.h"

namespace lupnt {

  // Mean and Osculating
  ClassicalOE Osculating2Mean(const ClassicalOE &coe_o, Real J2);
  ClassicalOE Mean2Osculating(const ClassicalOE &coe_m, Real J2);
  Vec6 Osculating2Mean(const Vec6 &coe_o, Real GM, Real J2);
  Vec6 Mean2Osculating(const Vec6 &coe_m, Real GM, Real J2);

  std::array<double, 6> ComputeSecondOrderShortPeriod(Vec6 &coe, Vec6 &doe);
  std::array<double, 6> ComputeFirstOrderMediumPeriod(Vec6 &coe, Vec6 &doe);
  std::array<double, 6> ComputeSecondOrderMediumPeriod(Vec6 &coe, Vec6 &doe);
  std::array<double, 6> ComputeCorrectionMediumPeriod(Vec6 &coe, Vec6 &doe);

}  // namespace lupnt
