/**
 * @file gnss_constellation.cpp
 * @author Stanford NAV LAB
 * @brief GNSS Constellation
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "gnss_constellation.h"

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

    double T = SECS_PER_DAY / tle.meanMotion;

    // Classical orbital elements
    ad::real a = pow((T * T * MU_EARTH) / (4.0 * PI * PI), 1.0 / 3.0);
    ad::real e = tle.eccentricity;
    ad::real i = tle.inclination * RAD_PER_DEG;
    ad::real Omega = tle.raan * RAD_PER_DEG;
    ad::real w = tle.argPerigee * RAD_PER_DEG;
    ad::real rad_per_sec = tle.meanMotion * 2 * PI /
                           SECS_PER_DAY;  // TLE mean motion is in revs/day
    ad::real M =
        tle.meanAnomaly * RAD_PER_DEG + dt_epoch * rad_per_sec;  // in radians

    ClassicalOE coe(a, e, i, Omega, w, M);
    coe.SetCoordSystem(CoordSystem::GCRF);
    auto cartOrbitState =
        std::make_shared<CartesianOrbitState>(CoeToCart(coe, MU_EARTH));

    // Create the spacecraft
    auto sat = std::make_shared<Spacecraft>();
    sat->SetDynamics(dynamics_);
    sat->SetOrbitState(cartOrbitState);
    sat->SetEpoch(epoch_);
    sat->SetBodyId(BodyId::EARTH);

    if (channel_) {
      auto transmitter = std::make_shared<GNSSTransmitter>(tle.name, tle.prn);
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
