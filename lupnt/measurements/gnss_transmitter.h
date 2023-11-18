/**
 * @file gnss_transmitter.h
 * @author Stanford NAV LAB
 * @brief  Handles Gnss transmit
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <lupnt/agents/comm_device.h>
#include <lupnt/measurements/transmission.h>

#include <string>
#include <tuple>

#include "antenna.h"
#include "gnss_measurement.h"

namespace lupnt {

class GnssChannel;

class GnssTransmitter : public Transmitter {
 public:
  Antenna antenna_;                    // Antenna gain pattern [deg & dB]
  std::string gnss_type_;              // Name of the atenna system
  std::string txrx = "TX";             // Type of comms system
  int prn_;                            // PRN of the transmitter satellite
  double freq_tx;                      // Transmit frequency [Hz]
  std::vector<std::string> freq_list;  // List of frequencies (by signal names)
  std::map<std::string, double> freq_map;  // map string to frequencies
  // Tramsmitter

  GnssTransmitter(std::string gnss_type, int prn)
      : gnss_type_(gnss_type), prn_(prn) {
    freq_map = {{"L1", 1575.42e6},  {"L2", 1227.60e6}, {"L5", 1176.45e6},
                {"E1", 1575.42e6},  {"E6", 1278.75e6}, {"E5", 1191.795e6},
                {"E5a", 1176.45e6}, {"E5b", 1207.14e6}};  // frequency maps
    InitializeGnssTransmitter();
  }

  // Transmitter Initializations
  void InitializeGnssTransmitter();
  void InitializeGPSTransmitter();
  void InitializeGLONASSTransmitter();
  void InitializeGALILEOTransmitter();
  void InitializeBEIDOUTransmitter();

  // Get transmitter orientatiion
  std::vector<Eigen::Vector3d> GetTransmitterOrientation(
      double t, Eigen::Vector3d& rv_tx_gcrf);
  double GetTransmittionAntennaGain(double t, Eigen::Vector3d r_tx_gcrf,
                                    Eigen::Vector3d r_rx_gcrf);

  // Get the Transmittion Information
  Transmission GenerateTransmission(double t);

  // Getters and Setters
  void SetChannel(std::shared_ptr<GnssChannel> ch) { channel = ch; };
  int GetPRN() { return prn_; };
  void SetFreq(double freq) { freq_tx = freq; };
  std::string GetGnssType() { return gnss_type_; };
  inline std::shared_ptr<Agent> GetAgent() const override { return agent; };
  inline void SetAgent(const std::shared_ptr<Agent>& agent) override {
    this->agent = agent;
  };
  double GetAntennaGain(Eigen::Vector3d direction) {
    return antenna_.GetAntennaGain(direction);
  };
  double GetAntennaGain(double theta, double phi) {
    return antenna_.GetAntennaGain(theta, phi);
  };

 private:
  std::shared_ptr<Agent> agent;  // Agent that owns the device
  std::shared_ptr<GnssChannel>
      channel;  // Channel that the device is connected to
};
}  // namespace lupnt
