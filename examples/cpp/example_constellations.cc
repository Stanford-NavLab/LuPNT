/**
 * @file ExampleConstellations.cpp
 * @author Stanford NAVLAB
 * @brief Setup Gnss Constellation
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
#include <lupnt/measurements/transmission.h>
#include <lupnt/numerics/filters.h>
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

int main() {
  auto cartesian_two_body =
      std::make_shared<CartesianTwoBodyDynamics>(GM_EARTH);

  auto channel = std::make_shared<GnssChannel>();

  // Gnss constellation
  for (std::string const_name : {"gps", "galileo", "beidou", "glonass"}) {
    auto gnss_const = GnssConstellation();
    gnss_const.SetDynamics(cartesian_two_body);
    gnss_const.SetChannel(channel);
    gnss_const.LoadTleFile(const_name);

    // Time
    double epoch_tai = gnss_const.GetEpoch();
    double t = epoch_tai;
    double dt = 1.0 * SECS_PER_MINUTE;
    double tf = t + 24.0 * SECS_PER_HOUR;
    double printEvery = 1.0 * SECS_PER_HOUR;

    // Output
    DataHistory data_history;
    FileWriter writer(std::filesystem::current_path() / "output" /
                          "ExampleConstellations" / const_name,
                      true);

    // Main loop
    while (t < tf) {
      gnss_const.Propagate(t);

      // GPS constellation
      for (int i = 0; i < gnss_const.GetNumSatellites(); i++) {
        auto sat_name = "sat" + std::to_string(i);
        auto sat_state = gnss_const.GetSatellite(i)
                             ->GetOrbitState()
                             ->GetVec()
                             .head(3)
                             .cast<double>();
        data_history.AddData(sat_name, t, sat_state);
      }

      // Print progress
      if (fmod(t, printEvery) < 1e-3) {
        std::cout << "t = " << std::fixed << std::setprecision(1) << t / 3600.0
                  << "/" << tf / 3600.0 << " hours" << std::endl;
      }

      t += dt;
    }

    // Write data
    writer.WriteData(data_history);
    std::cout << "Finished " << const_name << std::endl;
  }
  return 0;
}
