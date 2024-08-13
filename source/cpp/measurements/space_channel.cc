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

#include "lupnt/measurements/space_channel.h"

#include <iomanip>
#include <iostream>
#include <string>

#include "lupnt/agents/agent.h"
#include "lupnt/agents/comm_device.h"

namespace lupnt {

  void SpaceChannel::computeLinkBudget(std::shared_ptr<ICommDevice> &txDevice,
                                       std::shared_ptr<ICommDevice> &rxDevice, double t_tx,
                                       double t_rx, Transmission &transmission) {
    return;
  }
}  // namespace lupnt
