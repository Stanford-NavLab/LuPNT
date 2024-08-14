/**
 * @file gnss_channel.h
 * @author Stanford NAV LAB
 * @brief Gnss Channel class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <vector>

#include "lupnt/measurements/occultation.h"
#include "lupnt/measurements/space_channel.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {

  class GnssTransmitter;
  class GnssReceiver;

  class GnssChannel : public SpaceChannel {
  public:
    void AddTransmitter(std::shared_ptr<GnssTransmitter> &dev) { tx_devices.push_back(dev); }

    void AddReceiver(std::shared_ptr<GnssReceiver> &dev) { rx_devices.push_back(dev); }

    std::vector<Transmission> Receive(GnssReceiver &rx, double t);

    std::vector<std::shared_ptr<GnssReceiver>> rx_devices;
    std::vector<std::shared_ptr<GnssTransmitter>> tx_devices;
  };
}  // namespace lupnt
