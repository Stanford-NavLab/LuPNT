/**
 * @file gnss_constellation.cpp
 * @author Stanford NAV LAB
 * @brief Gnss Constellation
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/agents/gnss_constellation.h"

#include "lupnt/physics/orbit_state/tle.h"

namespace lupnt {

  void GnssConstellation::LoadTleFile(std::string filename) {
    auto path = kTlePath / filename;
    path += ".txt";

    bool is_first = true;

    for (auto tle : TLE::FromFile(path)) {
      // extract epoch
      if (is_first) {
        epoch_ = tle.epochTAI;
      }
      double dt_epoch = epoch_ - tle.epochTAI;

      double T = SECS_DAY / tle.meanMotion;

      // Classical orbital elements
      Real a = pow((T * T * GM_EARTH) / (4.0 * PI * PI), 1.0 / 3.0);
      Real e = tle.eccentricity;
      Real i = tle.inclination * RAD;
      Real Omega = tle.raan * RAD;
      Real w = tle.argPerigee * RAD;
      Real rad_per_sec = tle.meanMotion * 2 * PI / SECS_DAY;    // TLE mean motion is in revs/day
      Real M = tle.meanAnomaly * RAD + dt_epoch * rad_per_sec;  // in radians

      ClassicalOE coe({a, e, i, Omega, w, M});
      coe.SetCoordSystem(Frame::GCRF);
      auto cartOrbitState = std::make_shared<CartesianOrbitState>(Classical2Cart(coe, GM_EARTH));

      // Create the spacecraft
      auto sat = std::make_shared<Spacecraft>();
      sat->SetDynamics(dynamics_);
      sat->SetOrbitState(cartOrbitState);
      sat->SetEpoch(epoch_);
      sat->SetBodyId(NaifId::EARTH);

      if (channel_) {
        auto transmitter = std::make_shared<GnssTransmitter>(tle.name, tle.prn);
        sat->AddDevice(transmitter);
        transmitter->SetAgent(sat);
        channel_->AddTransmitter(transmitter);
        transmitter->SetChannel(channel_);
      }

      satellites_.push_back(sat);

      is_first = false;
    }
  }

}  // namespace lupnt