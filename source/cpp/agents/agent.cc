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

  std::shared_ptr<CartesianOrbitState> Agent::GetCartesianGCRFStateAtEpoch(Real epoch,
                                                                           Frame frame) {
    auto state = std::make_shared<OrbitState>(*state_);
    if (epoch != epoch_) dynamics_->Propagate(*state, epoch_, epoch, 1.0 * SECS_MINUTE);
    // TODO
    Real GM = 0.0;  //  GetBodyData(bodyId_).GM;
    auto cartOrbitState = std::static_pointer_cast<CartesianOrbitState>(
        ConvertOrbitStateRepresentation(state, OrbitStateRepres::CARTESIAN, GM));
    return ConvertOrbitStateFrame(cartOrbitState, epoch, Frame::GCRF);
  }

};  // namespace lupnt
