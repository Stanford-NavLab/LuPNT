/**
 * @file utils.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2025-02-06
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "lupnt/environment/plasma/core/math_utils.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace pecsim {

  double switchon(double x, double a, double da) {
    // Function varies from 0 to 1 to within 0.1% as the arguement x passes from
    // (a-da) to (a+da). So "da" is the range over which the function turns on.
    double c = 3.4534 / da;
    return tanh(c * (x - a)) / 2.0 + 0.5;
  }

  double sign(double x) { return (x > 0) - (x < 0); }

  double sign(double a, double b) { return (b >= 0) ? std::abs(a) : -std::abs(a); }

  double log10_safe(double value) {
    if (value > 0) {
      std::runtime_error("log10: value must be greater than 0.");
    }
    return std::log10(value);
  }

  double amod(double x, double y) { return x - y * std::floor(x / y); }

  double distance(const std::vector<double>& p1, const std::vector<double>& p2) {
    if (p1.size() != p2.size()) {
      throw std::runtime_error("Vectors must be of the same size.");
    }
    double sum = 0.0;
    for (size_t i = 0; i < p1.size(); ++i) {
      sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    }
    return std::sqrt(sum);
  }

}  // namespace pecsim
