/**
 * @file gnss_receiver.h
 * @author Stanford NAV LAB
 * @brief Gnss receiver class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "gnss_channel.h"
#include "gnss_measurement.h"
#include "lupnt/agents/comm_device.h"
#include "lupnt/measurements/transmission.h"

namespace lupnt {

  class GnssChannel;

  class GnssReceiver : public Receiver {
  public:
    Antenna antenna_;
    GnssReceiverParam gnssr_param_;
    std::string receiver_name_;  // Name of the receiver system

    GnssReceiver(std::string receiver_name) : receiver_name_(receiver_name) {
      antenna_ = Antenna(receiver_name_);
      InitializeReceiverParams();
    };

    // Receiver Gain Calculators
    std::vector<Vec3d> GetReceiverOrientation(double t, Vec3d& r_rx_gcrf, std::string mode);
    double GetReceiverAntennaGain(double t, Vec3d r_tx_gcrf, Vec3d r_rx_gcrf, std::string mode);

    void InitializeReceiverParams();
    void SetCN0Threshold(double CN0threshold) { rx_param_.CN0threshold = CN0threshold; };

    // Generate Measurement
    GnssMeasurement GetMeasurement(double t);
    void SetChannel(std::shared_ptr<GnssChannel> ch) { channel = ch; };

    // Getters and Setters
    inline std::shared_ptr<Agent> GetAgent() const override { return agent; };
    inline void SetAgent(const std::shared_ptr<Agent>& agent) override { this->agent = agent; };
    double GetAntennaGain(Vec3d direction) { return antenna_.GetAntennaGain(direction); };
    double GetAntennaGain(double theta, double phi) { return antenna_.GetAntennaGain(theta, phi); };

  private:
    std::shared_ptr<Agent> agent;
    std::shared_ptr<GnssChannel> channel;
  };
}  // namespace lupnt
