/**
 * @file agent.cpp
 * @author Stanford NAV Lab
 * @brief Base class fo
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "agent.h"

namespace lupnt {

int Agent::id_counter_ = 0;

std::shared_ptr<CartesianOrbitState> Agent::GetCartesianGCRFStateAtEpoch(
    ad::real epoch, CoordSystem coord_sys) {
  auto state = state_->Clone();
  if (epoch != epoch_)
    dynamics_->Propagate(*state, epoch_, epoch, 1.0 * SECS_PER_MINUTE);
  double mu = GetBodyData(bodyId_).GM;
  auto cartOrbitState = std::static_pointer_cast<CartesianOrbitState>(
      ConvertOrbitStateRepresentation(state, OrbitStateRepres::CARTESIAN, mu));
  return ConvertOrbitStateCoordSystem(cartOrbitState, epoch, CoordSystem::GCRF);
}

};  // namespace lupnt