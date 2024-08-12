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
#include "lupnt/core/constants.h"
#include "lupnt/measurements/comm_device.h"

namespace lupnt {

void SpaceChannel::ComputeLinkBudget(std::shared_ptr<ICommDevice> &txDevice,
                                     std::shared_ptr<ICommDevice> &rxDevice,
                                     double t_tx, double t_rx,
                                     Transmission &transmission) {
  return;
}

double SpaceChannel::ComputeFreeSpaceLossdB(double dist, double lambda) {
  double path_loss = pow(4 * PI * dist / lambda, 2);
  double path_loss_dB = 10 * log10(path_loss);
  return path_loss_dB;
}

}  // namespace lupnt