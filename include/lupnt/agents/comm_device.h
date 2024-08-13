/**
 * @file comm_device.h
 * @author Stanford NAV Lab
 * @brief Communication devices
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <map>
#include <memory>

#include "lupnt/measurements/antenna.h"
#include "lupnt/measurements/transmission.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/orbit_state.h"

namespace lupnt {

  class Agent;

  class ICommDevice {
  public:
    virtual ~ICommDevice() = default;
    std::string txrx = "none";
    virtual inline std::shared_ptr<Agent> GetAgent() const = 0;
    virtual inline void SetAgent(const std::shared_ptr<Agent> &agent) = 0;
  };

  struct TransmitterParam {
    double P_tx;  // Transmit power [dBW]
  };

  class Transmitter : public ICommDevice {
  public:
    virtual ~Transmitter() = default;
    std::string txrx = "tx";
    TransmitterParam tx_param_;
    virtual inline std::shared_ptr<Agent> GetAgent() const = 0;
    virtual inline void SetAgent(const std::shared_ptr<Agent> &agent) = 0;
  };

  struct ReceiverParam {
    double Ts;            // System noise temp [K]
    double Ae;            // Attenuation due to atmosphere (should be negative) [dB]
    double Nf;            // Noise figure of receiver/LNA [dB]
    double L;             // Receiver implementation, A/D conversion losses [dB]
    double As;            // System losses, in front of LNA [dB]
    double CN0threshold;  // CN0 threshold for receiving signals [dB-Hz]
  };

  class Receiver : public ICommDevice {
  public:
    virtual ~Receiver() = default;
    std::string txrx = "rx";
    ReceiverParam rx_param_;
    virtual inline std::shared_ptr<Agent> GetAgent() const = 0;
    virtual inline void SetAgent(const std::shared_ptr<Agent> &agent) = 0;
  };

  class Tranceiver : public ICommDevice {
  public:
    virtual ~Tranceiver() = default;
    std::string txrx = "txrx";
    TransmitterParam tx_param_;
    ReceiverParam rx_param_;
    virtual inline std::shared_ptr<Agent> GetAgent() const = 0;
    virtual inline void SetAgent(const std::shared_ptr<Agent> &agent) = 0;
  };

}  // namespace lupnt
