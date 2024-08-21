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
#include "lupnt/measurements/space_channel.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/orbit_state.h"

namespace lupnt {

class Agent;

enum CarrierType {
  Residual,  // Residual carrier
  BPSK,      // Binary Phase Shift Keying
  QPSK,      // Quadrature Phase Shift Keying
  OQPSK,     // Offset Quadrature Phase Shift Keying
  GMSK,      // Gaussian Minimum Shift Keying
  GMSK_PN,   // GMSK with PN modulation
};

enum FrequencyBand { UHF, L, S, Cband, X, Ku, K, Ka };

class ICommDevice {
 public:
  virtual ~ICommDevice() = default;
  std::string txrx = "none";
  std::string name = "none";
  inline virtual std::shared_ptr<Agent> GetAgent() { return agent_; };
  inline virtual void SetAgent(const std::shared_ptr<Agent> &agent) {
    agent_ = agent;
  };
  inline virtual std::shared_ptr<SpaceChannel> GetChannel() {
    return channel_;
  };
  inline virtual void SetChannel(const std::shared_ptr<SpaceChannel> &channel) {
    channel_ = channel;
  };

 private:
  std::shared_ptr<Agent> agent_;
  std::shared_ptr<SpaceChannel> channel_;
};

struct TransmitterParam {
  double P_tx = 0.0;              // Transmit power [dBW]
  double turnaround_ratio = 1.0;  // Turnaround ratio
};

class Transmitter : public ICommDevice {
 public:
  Antenna antenna_;
  double freq_tx;    // Transmit frequency [Hz]
  double bandwidth;  // Bandwidth of the signal [Hz]
  virtual ~Transmitter() = default;
  std::string txrx = "tx";
  TransmitterParam tx_param_;
  virtual inline double GetTransmittionAntennaGain(double t, Vec3d r_tx_gcrf,
                                                   Vec3d r_rx_gcrf) = 0;
};

struct ReceiverParam {
  double Tsys = 290;  // System noise temp [K]
  double Ae = -0.0;   // Attenuation due to atmosphere (should be negative) [dB]
  double L = -0.0;  // Receiver implementation, A/D conversion losses (should be
                    // negative) [dB]
  double As = -0.0;  // System losses, in front of LNA (should be negative) [dB]
  double CN0threshold = 20.0;  // CN0 threshold for receiving signals [dB-Hz]

  // PN Code Parameters
  CarrierType carrier_type = CarrierType::BPSK;  // carrier type
  double B_L_chip = 0.1;      // tracking loop noise bandwidth [Hz]
  double Tc = 1 / 2.068e6;    // chip duration
  double B_L_carrier = 0.1;   // carrier loop noise bandwidth [Hz]
  double m_R = 0.0;           // modulation index
  double T_I_doppler = 10.0;  // Doppler integration time [s]
  double T_I_range = 0.5;     // range integration time [s] (for open loop)
};

class Receiver : public ICommDevice {
 public:
  Antenna antenna_;
  virtual ~Receiver() = default;
  std::string txrx = "rx";
  ReceiverParam rx_param_;
  virtual inline double GetReceiverAntennaGain(double t, Vec3d r_tx_gcrf,
                                               Vec3d r_rx_gcrf) = 0;
};

class Transceiver : public ICommDevice {
 public:
  Antenna antenna_;
  double freq_tx;  // Transmit frequency [Hz]
  virtual ~Transceiver() = default;
  std::string txrx = "txrx";
  inline void SetTransmitter(const std::shared_ptr<Transmitter> &tx) {
    tx_ = tx;
  };
  inline void SetReceiver(const std::shared_ptr<Receiver> &rx) { rx_ = rx; };
  inline void SetAgent(const std::shared_ptr<Agent> &agent) override {
    SetAgent(agent);
    tx_->SetAgent(agent);
    rx_->SetAgent(agent);
  };

  inline std::shared_ptr<Transmitter> GetTransmitter() { return tx_; };
  inline std::shared_ptr<Receiver> GetReceiver() { return rx_; };

 private:
  std::shared_ptr<Transmitter> tx_;
  std::shared_ptr<Receiver> rx_;
};

}  // namespace lupnt
