#pragma once

#include "lupnt/states/state.h"

namespace lupnt {

  // Mean and Osculating
  ClassicalOE OsculatingToMean(const ClassicalOE &coe_o, Real GM, Real J2);
  ClassicalOE MeanToOsculating(const ClassicalOE &coe_m, Real GM, Real J2);
  Vec6 OsculatingToMean(const Vec6 &coe_o, Real GM, Real J2);
  Vec6 MeanToOsculating(const Vec6 &coe_m, Real GM, Real J2);

}  // namespace lupnt
