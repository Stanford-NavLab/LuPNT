/**
 * @file test_neldermead.cpp
 * @author Keidai Iiyama
 * @brief Test the Nelder-Mead optimization algorithm in 2D
 * @version 0.1
 * @date 2025-02-17
 */

#include <cmath>
#include <iostream>

#include "lupnt/environment/plasma/plasma.h"

using namespace pecsim;

int main() {
  // The ray-correction code can use Nelder-Mead; this toy Rosenbrock problem
  // provides a quick optimizer smoke test independent of the plasma model.
  auto rosenbrock = [](Vec2d p) {
    double x = p[0], y = p[1];
    return std::pow(1 - x, 2) + 100 * std::pow(y - x * x, 2);
  };

  NelderMead nm(rosenbrock, {0.0, 0.0});
  const Vec2d minimum = nm.minimize(2000, 1e-10, 1e-12, false);
  std::cout << "Minimum at: (" << minimum[0] << ", " << minimum[1] << "), f=" << nm.get_min_val()
            << "\n";
}
