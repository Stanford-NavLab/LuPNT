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

#include "lupnt/measurements/transmission.h"
#include "lupnt/physics/frame_converter.h"
#include "lupnt/physics/occultation.h"

namespace lupnt {

  struct Transmission;
  class GnssTransmitter;
  class GnssReceiver;

  class GnssChannel {
  public:
    void AddTransmitter(Ptr<GnssTransmitter> &dev) { tx_devices.push_back(dev); }

    void AddReceiver(Ptr<GnssReceiver> &dev) { rx_devices.push_back(dev); }

    std::vector<GnssTransmission> Receive(GnssReceiver &rx, double t);

    std::vector<Ptr<GnssReceiver>> rx_devices;
    std::vector<Ptr<GnssTransmitter>> tx_devices;
  };
}  // namespace lupnt
