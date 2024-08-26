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
#include "lupnt/measurements/comm_utils.h"
#include "lupnt/measurements/space_channel.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/orbit_state.h"

namespace lupnt {

  class Agent;

  class ICommDevice {
  public:
    virtual ~ICommDevice() = default;
    std::string txrx = "none";
    std::string name = "none";
    inline virtual void SetAgent(const std::shared_ptr<Agent> &agent) { agent_ = agent; };
    inline virtual Ptr<Agent> GetAgent() const { return agent_; };
    inline virtual Ptr<SpaceChannel> GetChannel() { return channel_; };
    inline virtual void SetChannel(const Ptr<SpaceChannel> &channel) { channel_ = channel; };

  protected:
    Ptr<Agent> agent_;
    Ptr<SpaceChannel> channel_;
  };

  class Transmitter : public ICommDevice {
  public:
    Antenna antenna_;
    double P_tx;       // Transmit power [dBW]
    double freq_tx;    // Transmit frequency [Hz]
    double bandwidth;  // Bandwidth of the signal [Hz]
    std::string txrx = "tx";
    Vec3 antenna_orientation_body = Vec3::Zero();

    virtual double GetTransmitterAntennaGain(double t, Vec3d r_tx_gcrf, Vec3d r_rx_gcrf) = 0;

    inline void SetAntennaOrientation(Vec3 orientation) { antenna_orientation_body = orientation; };
  };

  struct ReceiverParam {
    double Tsys = 290;           // System noise temp [K]
    double Ae = -0.0;            // Attenuation due to atmosphere (should be negative) [dB]
    double L = -0.0;             // Receiver implementation, A/D conversion losses (should be
                                 // negative) [dB]
    double As = -0.0;            // System losses, in front of LNA (should be negative) [dB]
    double CN0threshold = 20.0;  // CN0 threshold for receiving signals [dB-Hz]

    // PN Code Parameters
    Modulation modulation_type = Modulation::BPSK;  // carrier type
    double B_L_chip = 0.1;                          // tracking loop noise bandwidth [Hz]
    double Tc = 1 / 2.068e6;                        // chip duration
    double B_L_carrier = 0.1;                       // carrier loop noise bandwidth [Hz]
    double m_R = 0.0;                               // modulation index
    double T_I_doppler = 10.0;                      // Doppler integration time [s]
    double T_I_range = 0.5;                         // range integration time [s] (for open loop)
    std::string pn_ranging_code = "none";           //  "T2B", "T4B"
    double SER_threshold = 0.1;                     // Symbol error rate threshold
    double BTs = 0.5;            // For GMSK modulation (B: 3dB point of gaussian filter, )
    double coding_rate = 1 / 2;  // Coding rate
  };

  class Receiver : public ICommDevice {
  public:
    Antenna antenna_;
    std::string txrx = "rx";
    ReceiverParam rx_param_;
    Vec3 antenna_orientation_body = Vec3::Zero();
    virtual double GetReceiverAntennaGain(double t, Vec3d r_tx_gcrf, Vec3d r_rx_gcrf) = 0;

    inline void SetAntennaOrientation(Vec3 orientation) { antenna_orientation_body = orientation; };
  };

  class Transponder : public ICommDevice {
  public:
    Antenna antenna_;

    // transmitte params
    double freq_tx;                 // Transmit frequency [Hz]
    double turnaround_ratio = 1.0;  // Turnaround ratio
    std::string txrx = "txrx";

    Vec3 antenna_orientation_body = Vec3::Zero();

    Transponder(const std::shared_ptr<Transmitter> &tx, const std::shared_ptr<Receiver> &rx) {
      tx_ = tx;
      rx_ = rx;
    };
    virtual ~Transponder() = default;
    inline void SetTransmitter(const std::shared_ptr<Transmitter> &tx) { tx_ = tx; };
    inline void SetReceiver(const std::shared_ptr<Receiver> &rx) { rx_ = rx; };
    inline void SetAgent(const std::shared_ptr<Agent> &agent) override {
      SetAgent(agent);
      tx_->SetAgent(agent);
      rx_->SetAgent(agent);
    };
    inline void SetAntennaOrientation(Vec3 orientation) {
      antenna_orientation_body = orientation;
      tx_->SetAntennaOrientation(orientation);
      rx_->SetAntennaOrientation(orientation);
    };

    inline std::shared_ptr<Transmitter> GetTransmitter() { return tx_; };
    inline std::shared_ptr<Receiver> GetReceiver() { return rx_; };

  private:
    std::shared_ptr<Transmitter> tx_;
    std::shared_ptr<Receiver> rx_;
  };

}  // namespace lupnt
