#pragma once

#include <lupnt/lupnt.h>

namespace filtering_sim {
  using namespace lupnt;

  Real ComputeShapiroDelay(double t_tai_rx, const Vec3& tx_pos, const Vec3& rx_pos, Frame frame);
  Real ComputeRelativisticDelayLT(const Vec6& rv_sat, Frame frame, double dt);

}  // namespace filtering_sim
