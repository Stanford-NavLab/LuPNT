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
#include <lupnt/agents/agent.h>
#include <lupnt/agents/gnss_constellation.h>
#include <lupnt/agents/state_estimation_app_gnss.h>
#include <lupnt/core/file.h>
#include <lupnt/core/scheduler.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/measurements/gnss_channel.h>
#include <lupnt/measurements/gnss_receiver.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/clock.h>
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

int main() {
  // Time
  double t0 = 0;
  double tf = t0 + 36.0 * SECS_HOUR;  // 36 hours
  double dt = 1.0;                    // Integration time step [s]
  double Dt = 10.0;                   // Propagation time step [s]
  double print_every = 1.0 * SECS_HOUR;
  double save_every = Dt;

  // Dynamics
  auto dyn_earth_tb = std::make_shared<CartesianTwoBodyDynamics>(GM_EARTH);
  auto dyn_moon_tb = std::make_shared<CartesianTwoBodyDynamics>(GM_MOON);
  auto dyn_moon_nb = std::make_shared<NBodyDynamics>();

  auto earth = Body::Earth();
  auto moon = Body::Moon(5, 5);
  dyn_moon_nb->AddBody(moon);
  dyn_moon_nb->AddBody(earth);
  dyn_moon_nb->SetPrimaryBody(moon);

  auto dyn_clk = ClockDynamics(ClockModel::kMicrosemiCsac);

  // GPS constellation
  auto channel = std::make_shared<GnssChannel>();
  auto gps_const = GnssConstellation();
  gps_const.SetChannel(channel);
  gps_const.SetDynamics(dyn_earth_tb);
  gps_const.LoadTleFile("gps");

  // Epoch
  double epoch0 = gps_const.GetEpoch();
  std::string epoch_string = TAItoStringUTC(epoch0, 3);
  std::cout << "Initial Epoch: " << epoch_string << std::endl;

  // Moon spacecraft
  Real a = 6541.4;
  Real e = 0.6;
  Real i = 65.5 * RAD;
  Real Omega = 0.0 * RAD;
  Real w = 90.0 * RAD;
  Real M = 0.0 * RAD;
  ClassicalOE coe_moon({a, e, i, Omega, w, M});
  coe_moon.SetCoordSystem(Frame::MOON_CI);

  auto cart_state_moon =
      std::make_shared<CartesianOrbitState>(Classical2Cart(coe_moon, GM_MOON));
  auto moon_sat = std::make_shared<Spacecraft>();
  auto receiver = std::make_shared<GnssReceiver>("moongpsr");

  moon_sat->AddDevice(receiver);
  moon_sat->SetDynamics(dyn_moon_nb);
  moon_sat->SetOrbitState(cart_state_moon);
  moon_sat->SetEpoch(epoch0);
  moon_sat->SetBodyId(NaifId::MOON);
  moon_sat->SetClockDynamics(dyn_clk);

  receiver->SetAgent(moon_sat);
  receiver->SetChannel(channel);
  channel->AddReceiver(receiver);

  // Output
  auto data_history = std::make_shared<DataHistory>();
  auto output_path =
      std::filesystem::current_path() / "output" / "ExampleEKFApp";
  FileWriter writer(output_path, true);

  // OrbitState estimation app
  auto est_app = GnssStateEstimationApp();
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
    auto sate_mi = ConvertOrbitStateFrame(state, epoch, Frame::MOON_CI);
    auto state_gcrf = ConvertOrbitStateFrame(state, epoch, Frame::GCRF);
    data_history->AddData("rv_moon_mi", t, sate_mi->GetVec());
    data_history->AddData("rv_moon_gcrf", t, state_gcrf->GetVec());

    // GPS constellation
    for (int i = 0; i < gps_const.GetNumSatellites(); i++) {
      auto sate =
          gps_const.GetSatellite(i)->GetCartesianGCRFStateAtEpoch(epoch);
      auto sate_mi = ConvertOrbitStateFrame(sate, epoch, Frame::MOON_CI);
      auto state_gcrf = ConvertOrbitStateFrame(sate, epoch, Frame::GCRF);

      std::string name = "sat" + std::to_string(i);
      data_history->AddData(name + "_mi", t, sate->GetVec());
      data_history->AddData(name + "_gcrf", t, state_gcrf->GetVec());
    }

    // Bodies
    Vec6 v6;
    v6.setZero();

    // Vec6 vec6_mi = ConvertFrame(epoch, VecX::Zero(6),
    //                                           Frame::GCRF, Frame::MOON_CI);
    // Vec6 vec6_gcrf = ConvertFrame(epoch, VecX::Zero(6),
    //                                             Frame::MOON_CI, Frame::GCRF);

    // data_history->AddData("earth_mi", t, vec6_mi);
    // data_history->AddData("moon_gcrf", t, vec6_gcrf);

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