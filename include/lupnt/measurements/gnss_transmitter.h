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

#include <string>
#include <tuple>

#include "antenna.h"
#include "gnss_measurement.h"
#include "lupnt/measurements/comm_device.h"
#include "lupnt/measurements/transmission.h"

namespace lupnt {

  class GnssChannel;

  class GnssTransmitter : public Transmitter {
  public:
    Antenna antenna_;                        // Antenna gain pattern [deg & dB]
    std::string gnss_type_;                  // Name of the atenna system
    std::string txrx = "TX";                 // Type of comms system
    int prn_;                                // PRN of the transmitter satellite
    double freq_tx;                          // Transmit frequency [Hz]
    double Rc;                               // Ranging chip rate [Hz]
    std::vector<std::string> freq_list;      // List of frequencies (by signal names)
    std::map<std::string, double> freq_map;  // map string to frequencies
    std::map<std::string, double> rc_map;    // map string to chip rates

    // Tramsmitter

    GnssTransmitter(std::string gnss_type, int prn) : gnss_type_(gnss_type), prn_(prn) {
      freq_map = {{"L1", 1575.42e6},  {"L2", 1227.60e6}, {"L5", 1176.45e6},
                  {"E1", 1575.42e6},  {"E6", 1278.75e6}, {"E5", 1191.795e6},
                  {"E5a", 1176.45e6}, {"E5b", 1207.14e6}};  // frequency maps

      // Information from https://gnss-sdr.org/docs/tutorials/gnss-signals/
      rc_map
          = {{"L1", 1.023e6}, {"L2", 0.5115e6}, {"L5", 10.23e6}, {"E1", 1.023e6}, {"E6", 0.5115e6},
             {"E5", 10.23e6}, {"E5a", 10.23e6}, {"E5b", 10.23e6}};  // chip rate maps
      InitializeGnssTransmitter();
    }

    // Transmitter Initializations
    void InitializeGnssTransmitter();
    void InitializeGPSTransmitter();
    void InitializeGLONASSTransmitter();
    void InitializeGALILEOTransmitter();
    void InitializeBEIDOUTransmitter();

    // Get transmitter orientatiion
    std::vector<Vec3d> GetTransmitterOrientation(double t, Vec3d& rv_tx_gcrf);

    double GetTransmitterAntennaGain(double t, Vec3d r_tx_gcrf, Vec3d r_rx_gcrf) override;

    // Get the Transmittion Information
    GnssTransmission GenerateTransmission(double t);

    // Getters and Setters
    void SetChannel(Ptr<GnssChannel> ch) { channel = ch; };
    int GetPRN() { return prn_; };
    void SetFreq(double freq) { freq_tx = freq; };
    std::string GetGnssType() { return gnss_type_; };
    inline Ptr<Agent> GetAgent() const override { return agent; };
    inline void SetAgent(const Ptr<Agent>& agent) override { this->agent = agent; };
    double GetAntennaGain(Vec3d direction) { return antenna_.GetAntennaGain(direction); };
    double GetAntennaGain(double theta, double phi) { return antenna_.GetAntennaGain(theta, phi); };

  private:
    Ptr<Agent> agent;          // Agent that owns the device
    Ptr<GnssChannel> channel;  // Channel that the device is connected to
  };
}  // namespace lupnt
