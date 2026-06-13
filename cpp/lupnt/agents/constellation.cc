#include "lupnt/agents/constellation.h"

#include <regex>

#include "lupnt/agents/satellite.h"
#include "lupnt/conversions/state_conversions.h"
#include "lupnt/devices/comms.h"
#include "lupnt/dynamics/numerical_orbit_dynamics.h"
#include "lupnt/states/state.h"
#include "lupnt/states/tle.h"

namespace lupnt {

  Constellation::Constellation() {}

  Constellation::Constellation(Config& config) {
    auto tle_name = config["tle"].as<std::string>();
    name_ = config["name"].as<std::string>();
    for (auto tle : TLE::FromFile(tle_name)) {
      // Classical Orbital Elements
      Real T = SECS_DAY / tle.mean_motion;
      Real a = pow((T * T * GM_EARTH) / (4.0 * PI * PI), 1.0 / 3.0);
      Real e = tle.eccentricity;
      Real i = tle.inclination * RAD;
      Real Omega = tle.raan * RAD;
      Real w = tle.arg_perigee * RAD;
      Real rad_per_sec = tle.mean_motion * 2 * PI / SECS_DAY;
      Real dt = GetLupntEpoch() - tle.epoch_tai;
      Real M = WrapToPi(tle.mean_anomaly * RAD + rad_per_sec * dt);

      // State and Dynamics
      State coe = ClassicalOE({a, e, i, Omega, w, M}, Frame::GCRF);
      State rv = ClassicalToCart(coe, GM_EARTH);
      auto dyn = MakePtr<CartesianTwoBodyDynamics>(GM_EARTH);
      Config attitude_config(config["attitude_dynamics"]);
      auto attitude_dyn = MakePtr<FixedPointingDynamics>(attitude_config);
      Real frequency = config["frequency"].as<Real>();
      std::string name = name_ + "/" + std::regex_replace(tle.name, std::regex(" "), "_");

      // Devices
      auto gnss_transmitter = MakePtr<Transmitter>();
      gnss_transmitter->SetName(name + "/gnss_transmitter");

      // Satellite
      auto sat = MakePtr<Satellite>();
      sat->SetTime(0.0);
      sat->SetState(rv);
      sat->SetAttitude(Attitude(Vec4(1.0, 0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0), Frame::GCRF));
      sat->SetName(name);
      sat->SetFrequency(frequency);
      sat->SetDynamics(std::move(dyn));
      sat->AddDevice(std::move(gnss_transmitter));
      sat->SetAttitudeDynamics(std::move(attitude_dyn));

      satellites_.push_back(std::move(sat));
    }
  }

  void Constellation::Setup() {
    for (auto& sat : satellites_) sat->Setup();
  }

  void Constellation::Step(Real time) {
    for (auto& sat : satellites_) sat->Step(time);
  }

  void Constellation::SetSimulation(Simulation* sim) {
    sim_ = sim;
    for (auto& sat : satellites_) sat->SetSimulation(sim);
  }

}  // namespace lupnt
