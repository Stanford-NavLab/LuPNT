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
#include <lupnt/physics/frame_converter.h>
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
  auto earthDynamics = std::make_shared<CartesianTwoBodyDynamics>(GM_EARTH);
  auto moonDynamics = std::make_shared<CartesianTwoBodyDynamics>(GM_MOON);

  double epoch = 0.0;
  auto channel = std::make_shared<GnssChannel>();

  // GPS constellation
  auto gpsConstellation = GnssConstellation();
  gpsConstellation.SetChannel(channel);
  gpsConstellation.SetDynamics(earthDynamics);
  gpsConstellation.SetEpoch(epoch);
  gpsConstellation.LoadTleFile("gps");

  // Moon spacecraft
  Real a = 6541.4;
  Real e = 0.6;
  Real i = 65.5 * RAD;
  Real Omega = 90 * RAD;
  Real w = 0.0 * RAD;
  Real M = 0.0 * RAD;
  ClassicalOE coeMoon({a, e, i, Omega, w, M});
  coeMoon.SetCoordSystem(Frame::MOON_CI);

  auto cartOrbitStateMoon =
      std::make_shared<CartesianOrbitState>(Classical2Cart(coeMoon, GM_MOON));
  std::shared_ptr<Spacecraft> moonSat1 = std::make_shared<Spacecraft>();
  std::shared_ptr<GnssReceiver> receiver =
      std::make_shared<GnssReceiver>("moongpsr");

  moonSat1->AddDevice(receiver);
  moonSat1->SetDynamics(moonDynamics);
  moonSat1->SetOrbitState(cartOrbitStateMoon);
  moonSat1->SetEpoch(epoch);
  moonSat1->SetBodyId(NaifId::MOON);

  receiver->SetAgent(moonSat1);
  receiver->SetChannel(channel);
  channel->AddReceiver(receiver);

  // Time
  double t = epoch;
  double dt = 1.0 * SECS_MINUTE;
  double tf = t + 36.0 * SECS_HOUR;
  double printEvery = 1.0 * SECS_HOUR;

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
    auto stateMi = ConvertOrbitStateFrame(state, t, Frame::MOON_CI);
    auto stateGcrf = ConvertOrbitStateFrame(state, t, Frame::GCRF);
    dataHistory.AddData("moonSatMi", t, stateMi->GetVec());
    dataHistory.AddData("moonSatGcrf", t, stateGcrf->GetVec());

    // GPS constellation
    for (int i = 0; i < gpsConstellation.GetNumSatellites(); i++) {
      auto sate =
          gpsConstellation.GetSatellite(i)->GetCartesianGCRFStateAtEpoch(t);
      auto stateMi = ConvertOrbitStateFrame(sate, t, Frame::MOON_CI);
      auto stateGcrf = ConvertOrbitStateFrame(sate, t, Frame::GCRF);

      std::string name = "sat" + std::to_string(i);
      dataHistory.AddData(name + "Mi", t, sate->GetVec());
      dataHistory.AddData(name + "Gcrf", t, stateGcrf->GetVec());
    }

    // Bodies
    Vec6 v6;
    v6.setZero();
    dataHistory.AddData("earthMi", t,
                        ConvertFrame(t, v6, Frame::GCRF, Frame::MOON_CI));
    dataHistory.AddData("moonGcrf", t,
                        ConvertFrame(t, v6, Frame::MOON_CI, Frame::GCRF));

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