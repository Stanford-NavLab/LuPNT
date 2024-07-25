/**
 * @file gnss_constellation.h
 * @author Stanford NAV LAB
 * @brief
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <memory>
#include <vector>

#include "agent.h"
#include "lupnt/dynamics/dynamics.h"
#include "lupnt/measurements/gnss_channel.h"
#include "lupnt/measurements/gnss_transmitter.h"
#include "lupnt/physics/orbit_state.h"

namespace lupnt {

  class GnssConstellation {
  private:
    std::vector<std::shared_ptr<Spacecraft>> satellites_;
    std::shared_ptr<NumericalOrbitDynamics> dynamics_;
    std::shared_ptr<GnssChannel> channel_;
    double epoch_;  // in TAI
    const std::filesystem::path kTlePath = GetDataPath() / "tle";

  public:
    // Setters
    void SetChannel(std::shared_ptr<GnssChannel> ch) { channel_ = ch; }
    void SetDynamics(std::shared_ptr<NumericalOrbitDynamics> dyn) { dynamics_ = dyn; }
    void SetEpoch(double ep) { epoch_ = ep; }
    double GetEpoch() { return epoch_; }

    // Getters
    int GetNumSatellites() { return satellites_.size(); }
    std::shared_ptr<Spacecraft> GetSatellite(int i) { return satellites_[i]; }
    std::shared_ptr<GnssChannel> GetChannel() { return channel_; }
    std::shared_ptr<NumericalOrbitDynamics> GetDynamics() { return dynamics_; }

    // Methods
    void Propagate(double epoch) {
      if (epoch == epoch_) return;
      for (auto sat : satellites_) sat->Propagate(epoch);
      epoch_ = epoch;
    }

    void LoadTleFile(std::string filename);
  };

}  // namespace lupnt
