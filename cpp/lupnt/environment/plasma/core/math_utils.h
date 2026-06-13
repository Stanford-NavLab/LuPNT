/**
 * @file core/math_utils.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2025-02-06
 *
 * @copyright Copyright (c) 2025
 *
 */

#pragma once

#include <cmath>
#include <functional>
#include <vector>

#include "lupnt/environment/plasma/gcpm/constants_gcpm.h"

namespace pecsim {

  /**
   * @brief  Switching function varies from 0 to 1 to within 0.1% as the arguement
   * x passes from (a-da) to (a+da).
   *
   * @param x  the independent variable
   * @param a  the center of the switch
   * @param da  the range over which the function turns on
   * @return double
   */
  double switchon(double x, double a, double da);

  /**
   * @brief Returns the sign of a number.
   * @param x  the number to check
   * @return double 1 if x is positive, -1 if x is negative, and 0 if x is zero.
   */
  double sign(double x);

  /**
   * @brief Returns the value of a with a sign of b.
   * If b is positive, returns the absolute value of a.
   * If b is negative, returns the negative absolute value of a.
   *
   * @param a  the value to return
   * @param b  the sign to apply
   * @return double  the value of a with the sign of b
   */
  double sign(double a, double b);

  /**
   * @brief Computes the base-10 logarithm of a value, ensuring the value is
   * greater than 0 to avoid domain errors.
   *
   * @param value  the value to compute the logarithm for
   * @return double  the base-10 logarithm of the value
   */
  double log10_safe(double value);

  /**
   * @brief Computes the modulus of x with respect to y.
   *
   * @param x  the value to compute the modulus for
   * @param y  the divisor
   * @return double  the result of x mod y
   */
  double amod(double x, double y);

  /**
   * @brief Computes the Euclidean distance between two points in n-dimensional
   * space.
   *
   * @param p1  the first point as a vector of doubles
   * @param p2  the second point as a vector of doubles
   * @return double  the Euclidean distance between p1 and p2
   * @throws std::runtime_error if the vectors are not of the same size
   */
  double distance(const std::vector<double>& p1, const std::vector<double>& p2);

}  // namespace pecsim
