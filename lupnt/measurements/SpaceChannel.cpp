/**
 * @file SpaceChannel.cpp
 * @author Stanford NAV LAB
 * @brief Base spacechannel class (Under devlopment )
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "SpaceChannel.h"

#include <lupnt/agents/Agent.h>
#include <lupnt/agents/CommDevice.h>

#include <autodiff/forward/real.hpp>
#include <iomanip>
#include <iostream>
#include <string>

namespace LPT {

void SpaceChannel::computeLinkBudget(std::shared_ptr<ICommDevice> &txDevice,
                                     std::shared_ptr<ICommDevice> &rxDevice,
                                     double t_tx, double t_rx,
                                     Transmission &transmission) {
  return;
}
}  // namespace LPT