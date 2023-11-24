/**
 * @file ExamplePropagation.cpp
 * @author Stanford NAV LAB
 * @brief Example propagating both GPS constellation and Lunar Satellite
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

// lupnt includes
#include <lupnt/agents/agent.h>
#include <lupnt/agents/gnss_constellation.h>
#include <lupnt/core/file.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/measurements/gnss_channel.h>
#include <lupnt/measurements/gnss_measurement.h>
#include <lupnt/measurements/gnss_receiver.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/coord_converter.h>
#include <lupnt/physics/orbit_state.h>
#include <lupnt/physics/spice_interface.h>

// Autodiff includes

// Eigen includes

#include <Eigen/QR>

// C++ includes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace lupnt;

// Define the data structure for the TLE information

int main() {
  // Dynamics
  auto earthDynamics = std::make_shared<CartesianTwoBodyDynamics>(MU_EARTH);
  auto moonDynamics = std::make_shared<CartesianTwoBodyDynamics>(MU_MOON);

  double epoch = 0.0;
  auto channel = std::make_shared<GnssChannel>();

  // GPS constellation
  auto gpsConstellation = GnssConstellation();
  gpsConstellation.SetChannel(channel);
  gpsConstellation.SetDynamics(earthDynamics);
  gpsConstellation.SetEpoch(epoch);
  gpsConstellation.LoadTleFile("gps");

  // Moon spacecraft
  real a = 6541.4;
  real e = 0.6;
  real i = 65.5 * RAD_PER_DEG;
  real Omega = 90 * RAD_PER_DEG;
  real w = 0.0 * RAD_PER_DEG;
  real M = 0.0 * RAD_PER_DEG;
  ClassicalOE coeMoon({a, e, i, Omega, w, M});
  coeMoon.SetCoordSystem(CoordSystem::MI);

  auto cartOrbitStateMoon =
      std::make_shared<CartesianOrbitState>(CoeToCart(coeMoon, MU_MOON));
  std::shared_ptr<Spacecraft> moonSat1 = std::make_shared<Spacecraft>();
  std::shared_ptr<GnssReceiver> receiver =
      std::make_shared<GnssReceiver>("moongpsr");

  moonSat1->AddDevice(receiver);
  moonSat1->SetDynamics(moonDynamics);
  moonSat1->SetOrbitState(cartOrbitStateMoon);
  moonSat1->SetEpoch(epoch);
  moonSat1->SetBodyId(BodyId::MOON);

  receiver->SetAgent(moonSat1);
  receiver->SetChannel(channel);
  channel->AddReceiver(receiver);

  // Time
  double t = epoch;
  double dt = 1.0 * SECS_PER_MINUTE;
  double tf = t + 36.0 * SECS_PER_HOUR;
  double printEvery = 1.0 * SECS_PER_HOUR;

  // Output
  DataHistory dataHistory;
  FileWriter writer(
      std::filesystem::current_path() / "output" / "ExamplePropagation", true);

  // Main loop
  while (t < tf) {
    moonSat1->Propagate(t);
    gpsConstellation.Propagate(t);

    // Pseudorange
    auto measall = receiver->GetMeasurement(t);
    auto meas = measall.ExtractSignal("L1");

    // Naviigation
    dataHistory.AddData("pseudorange", t, meas.GetPseudorange());
    dataHistory.AddData("CN0", t, meas.GetCN0());

    dataHistory.AddData("earthOccultation", t, meas.GetEarthOccultation());
    dataHistory.AddData("moonOccultation", t, meas.GetMoonOccultation());
    dataHistory.AddData("antennaOccultation", t, meas.GetMoonOccultation());
    dataHistory.AddData("ionosOccultation", t, meas.GetMoonOccultation());

    // Moon spacecraft
    auto state = moonSat1->GetCartesianGCRFStateAtEpoch(t);
    auto stateMi = ConvertOrbitStateCoordSystem(state, t, CoordSystem::MI);
    auto stateGcrf = ConvertOrbitStateCoordSystem(state, t, CoordSystem::GCRF);
    dataHistory.AddData("moonSatMi", t, stateMi->GetVector());
    dataHistory.AddData("moonSatGcrf", t, stateGcrf->GetVector());

    // GPS constellation
    for (int i = 0; i < gpsConstellation.GetNumSatellites(); i++) {
      auto sate =
          gpsConstellation.GetSatellite(i)->GetCartesianGCRFStateAtEpoch(t);
      auto stateMi = ConvertOrbitStateCoordSystem(sate, t, CoordSystem::MI);
      auto stateGcrf = ConvertOrbitStateCoordSystem(sate, t, CoordSystem::GCRF);

      std::string name = "sat" + std::to_string(i);
      dataHistory.AddData(name + "Mi", t, sate->GetVector());
      dataHistory.AddData(name + "Gcrf", t, stateGcrf->GetVector());
    }

    // Bodies
    dataHistory.AddData(
        "earthMi", t,
        CoordConverter::Convert(VectorXreal::Zero(6), t, CoordSystem::GCRF,
                                CoordSystem::MI));
    dataHistory.AddData(
        "moonGcrf", t,
        CoordConverter::Convert(VectorXreal::Zero(6), t, CoordSystem::MI,
                                CoordSystem::GCRF));

    // Print progress
    if (fmod(t, printEvery) < 1e-3) {
      std::cout << "t = " << std::fixed << std::setprecision(1) << t / 3600.0
                << "/" << tf / 3600.0 << " hours" << std::endl;
    }

    t += dt;
  }

  // Write data
  writer.WriteData(dataHistory);
  std::cout << "Simulation finished" << std::endl;
  return 0;
}
