/**
 * @file GNSSChannel.h
 * @author Stanford NAV LAB
 * @brief GNSS Channel class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <lupnt/measurements/Occultation.h>
#include <lupnt/measurements/SpaceChannel.h>
#include <lupnt/physics/CoordConverter.h>

#include <vector>

namespace LPT {

class GNSSTransmitter;
class GNSSReceiver;

class GNSSChannel : public SpaceChannel {
 public:
  void AddTransmitter(std::shared_ptr<GNSSTransmitter> &dev) {
    tx_devices.push_back(dev);
  }

  void AddReceiver(std::shared_ptr<GNSSReceiver> &dev) {
    rx_devices.push_back(dev);
  }

  std::vector<Transmission> Receive(GNSSReceiver &rx, double t);

  std::vector<std::shared_ptr<GNSSReceiver>> rx_devices;
  std::vector<std::shared_ptr<GNSSTransmitter>> tx_devices;
};
}  // namespace LPT