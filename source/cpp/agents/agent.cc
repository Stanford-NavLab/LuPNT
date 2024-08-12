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
#include "lupnt/agents/agent.h"

namespace lupnt {

int Agent::id_counter_ = 0;

std::shared_ptr<CartesianOrbitState> Spacecraft::GetCartesianGCRFStateAtEpoch(
    Real epoch) {
  std::shared_ptr<OrbitState> state = GetOrbitState();
  Real current_epoch = GetEpoch();
  std::shared_ptr<NumericalOrbitDynamics> dynamics =
      std::dynamic_pointer_cast<NumericalOrbitDynamics>(GetDynamics());

  if (epoch != current_epoch) {
    // set dt
    Real dt = (epoch - current_epoch) / 10;
    dynamics->Propagate(*state, current_epoch, epoch, dt);
  }
  // TODO
  Real GM = 0.0;  //  GetBodyData(bodyId_).GM;
  auto cartOrbitState = std::static_pointer_cast<CartesianOrbitState>(
      ConvertOrbitStateRepresentation(state, OrbitStateRepres::CARTESIAN, GM));
  return ConvertOrbitStateFrame(cartOrbitState, epoch, Frame::GCRF);
}

};  // namespace lupnt