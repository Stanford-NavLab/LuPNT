// lupnt includes
#include <lupnt/agents/Agent.h>
#include <lupnt/agents/GNSSConstellation.h>
#include <lupnt/core/File.h>
#include <lupnt/dynamics/Dynamics.h>
#include <lupnt/measurements/GNSSChannel.h>
#include <lupnt/measurements/GNSSMeasurement.h>
#include <lupnt/measurements/GNSSReceiver.h>
#include <lupnt/numerics/MathUtils.h>
#include <lupnt/numerics/eigenmvn.h>
#include <lupnt/physics/CoordConverter.h>
#include <lupnt/physics/OrbitState.h>
#include <lupnt/physics/SpiceInterface.h>

// Autodiff includes
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/QR>

// C++ includes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace LPT;

// Define the data structure for the TLE information

int main() {
  // Dynamics
  auto earthDynamics = std::make_shared<CartesianTwoBodyDynamics>(MU_EARTH);
  auto moonDynamics = std::make_shared<CartesianTwoBodyDynamics>(MU_MOON);

  double epoch = 0.0;
  auto channel = std::make_shared<GNSSChannel>();

  // GPS constellation
  auto gpsConstellation = GNSSConstellation();
  gpsConstellation.SetChannel(channel);
  gpsConstellation.SetDynamics(earthDynamics);
  gpsConstellation.SetEpoch(epoch);
  gpsConstellation.LoadTleFile("gps");

  // Moon spacecraft
  ad::real a = 6541.4;
  ad::real e = 0.6;
  ad::real i = 65.5 * RAD_PER_DEG;
  ad::real Omega = 90 * RAD_PER_DEG;
  ad::real w = 0.0 * RAD_PER_DEG;
  ad::real M = 0.0 * RAD_PER_DEG;
  ClassicalOE coeMoon(a, e, i, Omega, w, M);
  coeMoon.SetCoordSystem(CoordSystem::MI);

  auto cartOrbitStateMoon =
      std::make_shared<CartesianOrbitState>(CoeToCart(coeMoon, MU_MOON));
  std::shared_ptr<Spacecraft> moonSat1 = std::make_shared<Spacecraft>();
  std::shared_ptr<GNSSReceiver> receiver =
      std::make_shared<GNSSReceiver>("moongpsr");

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
  FileWriter writer(BASEPATH / "output" / "ExamplePropagation", true);

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
        CoordConverter::Convert(ad::VectorXreal::Zero(6), t, CoordSystem::GCRF,
                                CoordSystem::MI));
    dataHistory.AddData(
        "moonGcrf", t,
        CoordConverter::Convert(ad::VectorXreal::Zero(6), t, CoordSystem::MI,
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
