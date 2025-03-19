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

namespace lupnt {

  void IntegratorParams::CheckIntegratorParams() {
    if (max_iter < 1) {
      throw std::invalid_argument("max_iter must be greater than 0");
    }
    if (abstol < 0) {
      throw std::invalid_argument("abstol must be non-negative");
    }
    if (reltol < 0) {
      throw std::invalid_argument("reltol must be non-negative");
    }
  }

  /**
   * @brief One step of Runge-Kutta Integration
   *
   * @param f  The ODE function to propagate
   * @param t  Time
   * @param x  The state to propagate
   * @param dt  Timestep
   */
  VecX RK4::Step(const ODE& f, Real t, const VecX& x, Real dt) {
    // usleep(1000);
    // return x + x;

    /* Evaluate `f` (i.e., `dx`) at the 4 locations defined bx the RK4 method */
    VecX k_1 = f(t, x) * dt;

    Real t1 = t + dt / 2.0;
    VecX x1 = x + k_1 / 2.0;
    VecX k_2 = f(t1, x1) * dt;

    Real t2 = t + dt / 2.0;
    VecX x2 = x + k_2 / 2.0;
    VecX k_3 = f(t2, x2) * dt;

    Real t3 = t + dt;
    VecX x3 = x + k_3;
    VecX k_4 = f(t3, x3) * dt;

    /* Average the 4 derivatives to approtimate `dx` */
    VecX dx = (k_1 + k_2 * 2.0 + k_3 * 2.0 + k_4) / 6.0;

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
  VecX RK8::Step(const ODE& f, Real t, const VecX& x, Real dt) {
    // 1
    VecX k_1 = f(t, x) * dt;

    // 2
    Real t1 = t + dt * (4.0 / 27.0);
    VecX x1 = x + k_1 * (4.0 / 27.0);
    VecX k_2 = f(t1, x1) * dt;

    // 3
    Real t2 = t + dt * (2.0 / 9);
    VecX x2 = x + (1.0 / 18) * (k_1 + 3 * k_2);
    VecX k_3 = f(t2, x2) * dt;

    // 4
    Real t3 = t + dt * (1.0 / 3);
    VecX x3 = x + (1.0 / 12) * (k_1 + 3 * k_3);
    VecX k_4 = f(t3, x3) * dt;

    // 5
    Real t4 = t + dt * (1.0 / 2);
    VecX x4 = x + (1.0 / 8) * (k_1 + 3 * k_4);
    VecX k_5 = f(t4, x4) * dt;

    // 6
    Real t5 = t + dt * (2.0 / 3);
    VecX x5 = x + (1.0 / 54) * (13 * k_1 - 27 * k_3 + 42 * k_4 + 8 * k_5);
    VecX k_6 = f(t5, x5) * dt;

    // 7
    Real t6 = t + dt * (1.0 / 6);
    VecX x6 = x + (1.0 / 4320) * (389 * k_1 - 54 * k_3 + 966 * k_4 - 824 * k_5 + 243 * k_6);
    VecX k_7 = f(t6, x6) * dt;

    // 8
    Real t7 = t + dt;
    VecX x7
        = x + (1.0 / 20) * (-234 * k_1 + 81 * k_3 - 1164 * k_4 + 656 * k_5 - 122 * k_6 + 800 * k_7);
    VecX k_8 = f(t7, x7) * dt;

    // 9
    Real t8 = t + (5.0 / 6) * dt;
    VecX x8
        = x
          + (1.0 / 288)
                * (-127 * k_1 + 18 * k_3 - 678 * k_4 + 456 * k_5 - 9 * k_6 + 576 * k_7 + 4 * k_8);
    VecX k_9 = f(t8, x8) * dt;

    // 10
    Real t9 = t + dt;
    VecX x9 = x
              + (1.0 / 820)
                    * (1481 * k_1 - 81 * k_3 + 7104 * k_4 - 3376 * k_5 + 72 * k_6 - 5040 * k_7
                       - 60 * k_8 + 720 * k_9);
    VecX k_10 = f(t9, x9) * dt;

    /* Approtimate `dx` */
    VecX dx
        = (41 * k_1 + 27 * k_4 + 272 * k_5 + 27 * k_6 + 216 * k_7 + 216 * k_9 + 41 * k_10) / 840;

    return x + dx;
  }

  /**
   * @brief One step of Runge-Kutta-Fehlberg Integration
   *
   * @param f  The ODE function to propagate
   * @param t  Time
   * @param x  The state to propagate
   * @param dt  Timestep
   *
   * @ref
   * https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
   */
  VecX IRKF::Step(const ODE& f, Real t, const VecX& x, Real dt) {
    int n = x.size();

    VecX x_new_low(n), x_new_high(n);

    for (int iter = 0; iter < params_.max_iter; ++iter) {
      Update(f, t, x, dt, x_new_low, x_new_high);
      bool within_tolerance = ComputeRelError(x_new_low, x_new_high, dt);

      if (within_tolerance) break;
    }

    return x_new_low;
  }

  bool IRKF::ComputeRelError(const VecX& x_new_low, const VecX& x_new_high, Real dt) {
    // Compute the error norm
    Real err_norm = 0.0;
    double tol = 0.0;
    double error = 0.0;
    bool within_tolerance = true;
    double max_error = 0.0;

    // Non-conservative acceptance threshold
    double accept_thresh
        = std::pow(order_ + 1,
                   (order_ + 1) / order_);  // J.C. Butcher, Numerical Methods for
                                            // Ordinary Differential Equations, p291

    for (size_t i = 0; i < x_new_low.size(); ++i) {
      error = std::abs(x_new_high(i).val() - x_new_low(i).val());
      tol = std::max<double>(params_.reltol * std::abs(x_new_high(i).val()), params_.abstol);
      max_error = std::max(max_error, error / tol);
      if ((error / tol) > accept_thresh) {
        within_tolerance = false;
      }
    }

    // Update timestep
    double beta = 0.9;  // safety factor
    double s = beta * std::pow(1.0 / max_error, 1.0 / (order_ + 1));
    s = std::max(0.5, std::min(2.0, s));  // J.C. Butcher, Numerical Methods for
                                          // Ordinary Differential Equations, (371a)
    dt = s * dt;

    return within_tolerance;
  }

  /**
   * @brief One step of Runge-Kutta-Fehlberg 45 Integration
   *
   * @param f  The ODE function to propagate
   * @param t  Time
   * @param x  The state to propagate
   * @param dt  Timestep
   *
   * @ref
   * https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
   */

  void RKF45::Update(const ODE& f, Real t, const VecX& x, const Real dt, VecX& x_new_low,
                     VecX& x_new_high) {
    // Compute k1 to k6
    VecX k1 = f(t, x) * dt;
    VecX k2 = f(t + dt / 4.0, x + k1 * (1.0 / 4.0)) * dt;
    VecX k3 = f(t + dt * 3.0 / 8.0, x + k1 * (3.0 / 32.0) + k2 * (9.0 / 32.0)) * dt;
    VecX k4 = f(t + dt * 12.0 / 13.0,
                x + k1 * (1932.0 / 2197.0) - k2 * (7200.0 / 2197.0) + k3 * (7296.0 / 2197.0))
              * dt;
    VecX k5 = f(t + dt,
                x + k1 * (439.0 / 216.0) - k2 * 8.0 + k3 * (3680.0 / 513.0) - k4 * (845.0 / 4104.0))
              * dt;
    VecX k6 = f(t + dt / 2.0, x - k1 * (8.0 / 27.0) + k2 * 2.0 - k3 * (3544.0 / 2565.0)
                                  + k4 * (1859.0 / 4104.0) - k5 * (11.0 / 40.0))
              * dt;

    // 4th order solution
    x_new_low = x + k1 * (25.0 / 216.0) + k3 * (1408.0 / 2565.0) + k4 * (2197.0 / 4104.0)
                - k5 * (1.0 / 5.0);

    // 5th order solution
    x_new_high = x + k1 * (16.0 / 135.0) + k3 * (6656.0 / 12825.0) + k4 * (28561.0 / 56430.0)
                 - k5 * (9.0 / 50.0) + k6 * (2.0 / 55.0);
  }
}  // namespace lupnt
