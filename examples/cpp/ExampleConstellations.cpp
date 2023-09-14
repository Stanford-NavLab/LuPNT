// lupnt includes
#include <lupnt/agents/Agent.h>
#include <lupnt/agents/GNSSConstellation.h>
#include <lupnt/core/File.h>
#include <lupnt/dynamics/Dynamics.h>
#include <lupnt/measurements/GNSSChannel.h>
#include <lupnt/measurements/GNSSMeasurement.h>
#include <lupnt/measurements/GNSSReceiver.h>
#include <lupnt/measurements/Transmission.h>
#include <lupnt/numerics/Filters.h>
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
namespace sp = SpiceInterface;

int main() {
  auto cartesian_two_body =
      std::make_shared<CartesianTwoBodyDynamics>(MU_EARTH);

  auto channel = std::make_shared<GNSSChannel>();

  // GNSS constellation
  for (std::string const_name : {"gps", "galileo", "beidou", "glonass"}) {
    auto gnss_const = GNSSConstellation();
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
    FileWriter writer(
        BASEPATH / "output" / "ExampleConstellations" / const_name, true);

    // Main loop
    while (t < tf) {
      gnss_const.Propagate(t);

      // GPS constellation
      for (int i = 0; i < gnss_const.GetNumSatellites(); i++) {
        auto sat_name = "sat" + std::to_string(i);
        auto sat_state = toEigen(
            gnss_const.GetSatellite(i)->GetOrbitState()->GetVector().head(3));
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
