/**
 * @file Integrator.cpp
 * @author Stanford NAV LAB
 * @brief Integrator implementations
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/numerics/integrator.h"

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

namespace ad = autodiff;

namespace lupnt {
/**
 * @brief One step of Runge-Kutta Integration
 *
 * @param f  The ODE function to propagate
 * @param t  Time
 * @param x  The state to propagate
 * @param dt  Timestep
 */
ad::VectorXreal RK4::Step(const ODE &f, const ad::real t,
                          const ad::VectorXreal &x, const ad::real dt) {
  /* Evaluate `f` (i.e., `dx`) at the 4 locations defined bx the RK4 method */
  ad::VectorXreal k_1 = f(t, x) * dt;

  ad::real t1 = t + dt / 2.0;
  ad::VectorXreal x1 = x + k_1 / 2.0;
  ad::VectorXreal k_2 = f(t1, x1) * dt;

  ad::real t2 = t + dt / 2.0;
  ad::VectorXreal x2 = x + k_2 / 2.0;
  ad::VectorXreal k_3 = f(t2, x2) * dt;

  ad::real t3 = t + dt;
  ad::VectorXreal x3 = x + k_3;
  ad::VectorXreal k_4 = f(t3, x3) * dt;

  /* Average the 4 derivatives to approtimate `dx` */
  ad::VectorXreal dx = (k_1 + k_2 * 2.0 + k_3 * 2.0 + k_4) / 6.0;

  return x + dx;
}

/**
 * @brief One step of 8th order Runge-Kutta Integration
 *
 * @param f  The ODE function to propagate
 * @param t  Time
 * @param x  The state to propagate
 * @param dt  Timestep
 *
 * @ref
 * https://www.mathworks.com/matlabcentral/fileexchange/55431-runge-kutta-8th-order-integration
 */
ad::VectorXreal RK8::Step(const ODE &f, const ad::real t,
                          const ad::VectorXreal &x, const ad::real dt) {
  // 1
  ad::VectorXreal k_1 = f(t, x) * dt;

  // 2
  ad::real t1 = t + dt * (4.0 / 27.0);
  ad::VectorXreal x1 = x + k_1 * (4.0 / 27.0);
  ad::VectorXreal k_2 = f(t1, x1) * dt;

  // 3
  ad::real t2 = t + dt * (2.0 / 9);
  ad::VectorXreal x2 = x + (1.0 / 18) * (k_1 + 3 * k_2);
  ad::VectorXreal k_3 = f(t2, x2) * dt;

  // 4
  ad::real t3 = t + dt * (1.0 / 3);
  ad::VectorXreal x3 = x + (1.0 / 12) * (k_1 + 3 * k_3);
  ad::VectorXreal k_4 = f(t3, x3) * dt;

  // 5
  ad::real t4 = t + dt * (1.0 / 2);
  ad::VectorXreal x4 = x + (1.0 / 8) * (k_1 + 3 * k_4);
  ad::VectorXreal k_5 = f(t4, x4) * dt;

  // 6
  ad::real t5 = t + dt * (2.0 / 3);
  ad::VectorXreal x5 =
      x + (1.0 / 54) * (13 * k_1 - 27 * k_3 + 42 * k_4 + 8 * k_5);
  ad::VectorXreal k_6 = f(t5, x5) * dt;

  // 7
  ad::real t6 = t + dt * (1.0 / 6);
  ad::VectorXreal x6 = x + (1.0 / 4320) * (389 * k_1 - 54 * k_3 + 966 * k_4 -
                                           824 * k_5 + 243 * k_6);
  ad::VectorXreal k_7 = f(t6, x6) * dt;

  // 8
  ad::real t7 = t + dt;
  ad::VectorXreal x7 = x + (1.0 / 20) * (-234 * k_1 + 81 * k_3 - 1164 * k_4 +
                                         656 * k_5 - 122 * k_6 + 800 * k_7);
  ad::VectorXreal k_8 = f(t7, x7) * dt;

  // 9
  ad::real t8 = t + (5.0 / 6) * dt;
  ad::VectorXreal x8 =
      x + (1.0 / 288) * (-127 * k_1 + 18 * k_3 - 678 * k_4 + 456 * k_5 -
                         9 * k_6 + 576 * k_7 + 4 * k_8);
  ad::VectorXreal k_9 = f(t8, x8) * dt;

  // 10
  ad::real t9 = t + dt;
  ad::VectorXreal x9 =
      x + (1.0 / 820) * (1481 * k_1 - 81 * k_3 + 7104 * k_4 - 3376 * k_5 +
                         72 * k_6 - 5040 * k_7 - 60 * k_8 + 720 * k_9);
  ad::VectorXreal k_10 = f(t9, x9) * dt;

  /* Approtimate `dx` */
  ad::VectorXreal dx = (41 * k_1 + 27 * k_4 + 272 * k_5 + 27 * k_6 + 216 * k_7 +
                        216 * k_9 + 41 * k_10) /
                       840;

  return x + dx;
}

}  // namespace lupnt
