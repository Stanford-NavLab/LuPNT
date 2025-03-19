/**
 * @file ExamplePerformance.cpp
 * @author Stanford NAV LAB
 * @brief Example for performance testing of spherical harmonics computation
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

// lupnt includes
#include <lupnt/core/file.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/frame_converter.h>
#include <lupnt/physics/orbit_state.h>
#include <lupnt/physics/spice_interface.h>

// Autodiff includes

// Eigen includes

#include <Eigen/QR>

// C++ includes
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace lupnt;
using namespace std::chrono;

int main() {
  // Bodies
  Body earth = Body::Earth();

  for (int h : {0, 25, 50, 100}) {
    Body moon = Body::Moon();
    moon.gravity_field = ReadHarmonicGravityField<double>("grgm900c.cof", h, h, false);
    // moon.sphericalHarmonics = false;

    // Dynamics
    NBodyDynamics<double> dynamics;
    dynamics.AddBody(moon);
    // dynamics.AddBody(earth);
    // dynamics.SetPrimaryBody(moon);

    // State
    VecX rv0(6);
    rv0 << -1.540113643726188e3, -0.179443941906269e3, 1.128341549807345e3, -0.000291469032495e3,
        -0.001449961303523e3, -0.000628428693161e3;
    CartesianOrbitState cart_state(rv0);

    // Time
    double tStart = 662731269.183929;
    double dt = 10.0 * SECS_HOUR;
    double tEnd = tStart + 5.0 * SECS_DAY;
    double t = tStart;

    // NBodyRates
    int n = 100;
    int m = 24 * 60;
    VecXd times(n);
    VecX a;
    VecXd J;
    // std::function<VecX(real, VecX &, NBodyDynamics
    // &)>
    //     func = [=](real t, VecX &xs, NBodyDynamics &dyn) {
    //       return dyn.ComputeRates(t, xs);
    //     };

    for (int i = 0; i < n; i++) {
      auto start = high_resolution_clock::now();
      for (int j = 0; j < m; j++) {
        // J = jacobian(func, wrt(rv0), at(tStart, rv0, dynamics), a);
        dynamics.ComputeRates(tStart, rv0);
      }
      auto end = high_resolution_clock::now();
      auto duration = duration_cast<microseconds>(end - start);
      times(i) = duration.count();
      // std::cout << "Run: " << i << " Time: " << duration.count() / 1e6 << "
      // s"
      //           << std::endl;
    }
    double n_double = n;
    double mean = times.sum() / n_double;
    double std_dev = std::sqrt((times.array() - times.mean()).square().sum() / (n_double - 1.0));

    std::cout << "NBodyRates: " << h << std::endl;
    std::cout << "Mean, Std: " << mean / 1e6 << ", " << std_dev / 1e6 << " s" << std::endl;
    std::cout << "Total: " << times.sum() / 1e6 << " s" << std::endl;
    t += dt;
  }

  return 0;
}
