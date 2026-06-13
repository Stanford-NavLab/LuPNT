/**
 * @file ExamplePerformance.cpp
 * @author Stanford NAV LAB
 * @brief Benchmark spherical-harmonic gravity acceleration evaluation.
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <lupnt/lupnt.h>

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
  for (int h : {0, 25, 50, 100}) {
    // Use the body factory so the loaded gravity field matches Body's scalar
    // type and unit system.
    Body moon = Body::Moon(h, h);

    // Dynamics
    NBodyDynamics dynamics;
    dynamics.AddBody(moon);
    dynamics.SetFrame(Frame::MOON_CI);
    // dynamics.AddBody(earth);
    // dynamics.SetPrimaryBody(moon);

    // State
    Vec6 rv0{-1.540113643726188e3, -0.179443941906269e3, 1.128341549807345e3,
             -0.000291469032495e3, -0.001449961303523e3, -0.000628428693161e3};

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
