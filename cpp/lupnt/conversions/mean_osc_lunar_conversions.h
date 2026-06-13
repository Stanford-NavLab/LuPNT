#pragma once

#include "lupnt/states/state.h"

namespace lupnt {

  Vec6 LunarMeanToOsculating(Vec6 meanCoeVec);
  Vec6 LunarOsculatingToMean(Vec6 oscCoeVec);

  std::array<double, 6> ComputeSecondOrderShortPeriod(Vec6 &coe, Vec6 &doe);
  std::array<double, 6> ComputeFirstOrderMediumPeriod(Vec6 &coe, Vec6 &doe);
  std::array<double, 6> ComputeSecondOrderMediumPeriod(Vec6 &coe, Vec6 &doe);
  std::array<double, 6> ComputeCorrectionMediumPeriod(Vec6 &coe, Vec6 &doe);

}  // namespace lupnt
