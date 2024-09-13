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
#include "lupnt/measurements/comm_device.h"
#include "lupnt/measurements/transmission.h"

namespace lupnt {

  class GnssChannel;

  class GnssReceiver : public Receiver {
  public:
    Antenna antenna_;
    GnssReceiverParam gnssr_param_;
    std::string receiver_name_;           // Name of the receiver system
    std::string attitude_mode_ = "NONE";  // Attitude mode of the receiver

    GnssReceiver(std::string receiver_name) : receiver_name_(receiver_name) {
      antenna_ = Antenna(receiver_name_);
      InitializeReceiverParams();
    };

    // Receiver Gain Calculators
    void SetReceiverAttitudeMode(std::string mode) { attitude_mode_ = mode; };

    // Receiver Gain Calculators
    std::vector<Vec3d> GetReceiverOrientation(double t, Vec3d& r_rx_gcrf, std::string mode);

    /**
     * @brief Get the Receiver Antenna Gain object
     *
     * @param t   epoch (TAI) [s]
     * @param r_tx_gcrf  position of the transmitter in GCRF [km]
     * @param r_rx_gcrf  position of the receiver in GCRF [km]
     * @param frame_tx   frame of the transmitter
     * @param frame_rx   frame of the receiver
     *
     * @return double
     */
    double GetReceiverAntennaGain(double t, Vec3d r_tx_gcrf, Vec3d r_rx_gcrf) override;

    void InitializeReceiverParams();
    void SetCN0Threshold(double CN0threshold) { rx_param_.CN0threshold = CN0threshold; };

    // Generate Measurement
    GnssMeasurement GetMeasurement(double t);

    // Getters and Setters
    // Override channel getters and setters
    inline Ptr<Agent> GetAgent() const override { return agent; };
    inline void SetAgent(Ptr<Agent> agent) override { this->agent = agent; };
    inline Ptr<GnssChannel> GetGnssChannel() { return gnss_channel_; };
    inline void SetChannel(Ptr<SpaceChannel> channel) override {
      gnss_channel_ = std::static_pointer_cast<GnssChannel>(channel);
    };

    double ComputeGain(Vec3d direction) { return antenna_.ComputeGain(direction); };
    double ComputeGain(double theta, double phi) { return antenna_.ComputeGain(theta, phi); };

  private:
    Ptr<Agent> agent;  // Agent that owns the device
    Ptr<GnssChannel> gnss_channel_;
  };
}  // namespace lupnt
