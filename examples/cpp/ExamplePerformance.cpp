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
#include <lupnt/core/File.h>
#include <lupnt/dynamics/Dynamics.h>
#include <lupnt/dynamics/GravityField.h>
#include <lupnt/numerics/MathUtils.h>
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
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace LPT;
using namespace std::chrono;
namespace sp = SpiceInterface;

int main() {
  // Bodies
  Body earth = Body::Earth();

  for (int h : {0, 25, 50, 100}) {
    Body moon = Body::Moon(h, h);
    // moon.sphericalHarmonics = false;

    // Dynamics
    NBodyDynamics dynamics;
    dynamics.AddBody(moon);
    // dynamics.AddBody(earth);
    dynamics.SetCentralBody(moon);

    // State
    ad::VectorXreal rv0(6);
    rv0 << -1.540113643726188e3, -0.179443941906269e3, 1.128341549807345e3,
        -0.000291469032495e3, -0.001449961303523e3, -0.000628428693161e3;
    CartesianOrbitState cart_state(rv0);

    // Time
    double tStart = 662731269.183929;
    double dt = 10.0 * SECS_PER_HOUR;
    double tEnd = tStart + 5.0 * SECS_PER_DAY;
    double t = tStart;

    // NBodyRates
    int n = 100;
    int m = 24 * 60;
    Eigen::ArrayXd times(n);
    ad::VectorXreal a;
    Eigen::MatrixXd J;
    // std::function<ad::VectorXreal(ad::real, ad::VectorXreal &, NBodyDynamics
    // &)>
    //     func = [=](ad::real t, ad::VectorXreal &xs, NBodyDynamics &dyn) {
    //       return dyn.ComputeRates(t, xs);
    //     };

    for (int i = 0; i < n; i++) {
      auto start = high_resolution_clock::now();
      for (int j = 0; j < m; j++) {
        // J = ad::jacobian(func, wrt(rv0), at(tStart, rv0, dynamics), a);
        dynamics.ComputeRates(tStart, rv0);
      }
      auto end = high_resolution_clock::now();
      auto duration = duration_cast<microseconds>(end - start);
      times(i) = duration.count();
      // std::cout << "Run: " << i << " Time: " << duration.count() / 1e6 << "
      // s"
      //           << std::endl;
    }
    double mean = times.sum() / (double)times.size();
    double std_dev = std::sqrt((times - times.mean()).square().sum() /
                               (double)(times.size() - 1.0));

    std::cout << "NBodyRates: " << h << std::endl;
    std::cout << "Mean, Std: " << mean / 1e6 << ", " << std_dev / 1e6 << " s"
              << std::endl;
    std::cout << "Total: " << times.sum() / 1e6 << " s" << std::endl;
    t += dt;
  }

  return 0;
}
