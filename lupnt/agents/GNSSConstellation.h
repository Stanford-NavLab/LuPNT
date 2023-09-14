#pragma once

#include <lupnt/agents/Agent.h>
#include <lupnt/dynamics/Dynamics.h>
#include <lupnt/measurements/GNSSChannel.h>
#include <lupnt/measurements/GNSSTransmitter.h>
#include <lupnt/physics/OrbitState.h>

#include <memory>
#include <vector>
namespace LPT {

class GNSSConstellation {
 private:
  std::vector<std::shared_ptr<Spacecraft>> satellites_;
  std::shared_ptr<IOrbitDynamics> dynamics_;
  std::shared_ptr<GNSSChannel> channel_;
  double epoch_;  // in TAI
  const std::filesystem::path kTlePath = BASEPATH / "data" / "tle";

 public:
  // Setters
  void SetChannel(std::shared_ptr<GNSSChannel> ch) { channel_ = ch; }
  void SetDynamics(std::shared_ptr<IOrbitDynamics> dyn) { dynamics_ = dyn; }
  void SetEpoch(double ep) { epoch_ = ep; }
  double GetEpoch() { return epoch_; }

  // Getters
  int GetNumSatellites() { return satellites_.size(); }
  std::shared_ptr<Spacecraft> GetSatellite(int i) { return satellites_[i]; }
  std::shared_ptr<GNSSChannel> GetChannel() { return channel_; }
  std::shared_ptr<IOrbitDynamics> GetDynamics() { return dynamics_; }

  // Methods
  void Propagate(double epoch) {
    if (epoch == epoch_) return;
    for (auto sat : satellites_) sat->Propagate(epoch);
    epoch_ = epoch;
  }

  void LoadTleFile(std::string filename);
};

}  // namespace LPT
