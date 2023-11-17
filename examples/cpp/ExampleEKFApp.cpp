/**
 * @file ExampleEKFApp.cpp
 * @author Stanford NAV LAB
 * @brief Example of EKF application
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

// lupnt includes
#include <lupnt/agents/Agent.h>
#include <lupnt/agents/GNSSConstellation.h>
#include <lupnt/agents/GNSSStateEstimationApp.h>
#include <lupnt/core/File.h>
#include <lupnt/core/Scheduler.h>
#include <lupnt/dynamics/Dynamics.h>
#include <lupnt/measurements/GNSSChannel.h>
#include <lupnt/measurements/GNSSReceiver.h>
#include <lupnt/numerics/MathUtils.h>
#include <lupnt/physics/Clock.h>
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
namespace sp = SpiceInterface;

int main() {
  // Time
  double t0 = 0;
  double tf = t0 + 36.0 * SECS_PER_HOUR;  // 36 hours
  double dt = 1.0;                        // Integration time step [s]
  double Dt = 10.0;                       // Propagation time step [s]
  double print_every = 1.0 * SECS_PER_HOUR;
  double save_every = Dt;

  // Dynamics
  auto dyn_earth_tb = std::make_shared<CartesianTwoBodyDynamics>(MU_EARTH);
  auto dyn_moon_tb = std::make_shared<CartesianTwoBodyDynamics>(MU_MOON);
  auto dyn_moon_nb = std::make_shared<NBodyDynamics>();

  auto earth = Body::Earth();
  auto moon = Body::Moon(5, 5);
  dyn_moon_nb->AddBody(moon);
  dyn_moon_nb->AddBody(earth);
  dyn_moon_nb->SetCentralBody(moon);

  auto dyn_clk = ClockDynamics(ClockModel::kMicrosemiCsac);

  // GPS constellation
  auto channel = std::make_shared<GNSSChannel>();
  auto gps_const = GNSSConstellation();
  gps_const.SetChannel(channel);
  gps_const.SetDynamics(dyn_earth_tb);
  gps_const.LoadTleFile("gps");

  // Epoch
  double epoch0 = gps_const.GetEpoch();
  std::string epoch_string = sp::TAItoStringUTC(epoch0, 3);
  std::cout << "Initial Epoch: " << epoch_string << std::endl;

  // Moon spacecraft
  ad::real a = 6541.4;
  ad::real e = 0.6;
  ad::real i = 65.5 * RAD_PER_DEG;
  ad::real Omega = 0.0 * RAD_PER_DEG;
  ad::real w = 90.0 * RAD_PER_DEG;
  ad::real M = 0.0 * RAD_PER_DEG;
  ClassicalOE coe_moon(a, e, i, Omega, w, M);
  coe_moon.SetCoordSystem(CoordSystem::MI);

  auto cart_state_moon =
      std::make_shared<CartesianOrbitState>(CoeToCart(coe_moon, MU_MOON));
  auto moon_sat = std::make_shared<Spacecraft>();
  auto receiver = std::make_shared<GNSSReceiver>("moongpsr");

  moon_sat->AddDevice(receiver);
  moon_sat->SetDynamics(dyn_moon_nb);
  moon_sat->SetOrbitState(cart_state_moon);
  moon_sat->SetEpoch(epoch0);
  moon_sat->SetBodyId(BodyId::MOON);
  moon_sat->SetClockDynamics(dyn_clk);

  receiver->SetAgent(moon_sat);
  receiver->SetChannel(channel);
  channel->AddReceiver(receiver);

  // Output
  auto data_history = std::make_shared<DataHistory>();
  auto output_path = BASEPATH / "output" / "ExampleEKFApp";
  FileWriter writer(output_path, true);

  // OrbitState estimation app
  auto est_app = GNSSStateEstimationApp();
  est_app.SetAgent(moon_sat);
  est_app.SetDynamics(dyn_moon_tb);
  est_app.SetReceiver(receiver);
  est_app.SetDataHistory(data_history);
  est_app.SetFrequency(10.0);

  // Propagate and save function
  auto PropagateAndSave = [&](double t) {
    double epoch = epoch0 + t;

    // Propagate
    moon_sat->Propagate(epoch);
    gps_const.Propagate(epoch);

    // Moon spacecraft
    auto state = moon_sat->GetCartesianGCRFStateAtEpoch(epoch);
    auto sate_mi = ConvertOrbitStateCoordSystem(state, epoch, CoordSystem::MI);
    auto state_gcrf =
        ConvertOrbitStateCoordSystem(state, epoch, CoordSystem::GCRF);
    data_history->AddData("rv_moon_mi", t, sate_mi->GetVector());
    data_history->AddData("rv_moon_gcrf", t, state_gcrf->GetVector());

    // GPS constellation
    for (int i = 0; i < gps_const.GetNumSatellites(); i++) {
      auto sate =
          gps_const.GetSatellite(i)->GetCartesianGCRFStateAtEpoch(epoch);
      auto sate_mi = ConvertOrbitStateCoordSystem(sate, epoch, CoordSystem::MI);
      auto state_gcrf =
          ConvertOrbitStateCoordSystem(sate, epoch, CoordSystem::GCRF);

      std::string name = "sat" + std::to_string(i);
      data_history->AddData(name + "_mi", t, sate->GetVector());
      data_history->AddData(name + "_gcrf", t, state_gcrf->GetVector());
    }

    // Bodies
    data_history->AddData(
        "earth_mi", t,
        CoordConverter::Convert(ad::VectorXreal::Zero(6), epoch,
                                CoordSystem::GCRF, CoordSystem::MI));
    data_history->AddData(
        "moon_gcrf", t,
        CoordConverter::Convert(ad::VectorXreal::Zero(6), epoch,
                                CoordSystem::MI, CoordSystem::GCRF));

    // Print progress
    if (fmod(t, print_every) < 1e-3) {
      std::cout << "t = " << std::fixed << std::setprecision(1) << t / 3600.0
                << "/" << tf / 3600.0 << " hours (" << (t / tf * 100) << "%)"
                << std::endl
                << std::setprecision(4);
    }
  };

  // Start time
  auto start = std::chrono::high_resolution_clock::now();

  // Scheduler
  Scheduler scheduler;
  scheduler.ScheduleFunction(PropagateAndSave, t0 + Dt,
                             Dt);  // frequency with Dt
  scheduler.ScheduleApplication(est_app, t0 + Dt, Dt);
  scheduler.RunSimulation(tf);

  // double t = t0;
  // est_app.Setup();
  // while (t < tf) {
  //   t += Dt;
  //   PropagateAndSave(t);
  //   est_app.Step(t);
  // }

  // End time
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;

  // Write data
  std::cout << "Simulation finished in " << elapsed_seconds.count()
            << " s, saving data..." << std::endl;
  writer.WriteData(*data_history);
  std::cout << "Data saved to " << output_path << std::endl;
  return 0;
}
