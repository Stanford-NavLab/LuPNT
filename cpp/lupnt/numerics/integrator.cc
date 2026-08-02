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

#include "lupnt/core/error.h"
#include "lupnt/core/progress_bar.h"

namespace lupnt {

  void IntegratorParams::CheckIntegratorParams() {
    LUPNT_CHECK(max_iter > 0, "max_iter must be greater than 0", "IntegratorParams");
    LUPNT_CHECK(abstol > 0, "abstol must be non-negative", "IntegratorParams");
    LUPNT_CHECK(reltol > 0, "reltol must be non-negative", "IntegratorParams");
  }

  State Integrator::Propagate(const ODE& odefunc, Real t0, Real tf, const State& x0, Real dt) {
    LUPNT_CHECK(dt > 0, "dt must be positive", "Integrator");
    State x = x0;
    Real t = t0;

    Ptr<ProgressBar> pbar = nullptr;
    if (print_progress_) pbar = Logger::GetProgressBar(100, "", "Integrator");
    while (t < tf) {
      dt = std::min(dt, tf - t);
      Real prev_dt = dt;
      x = Step(odefunc, t, x, dt);
      t += prev_dt;
      if (print_progress_) pbar->Update(int((t - t0) / (tf - t0) * 100));
    }
    if (print_progress_) pbar->Finish();
    return x;
  };

  State Integrator::Propagate(const ODE& odefunc, Real t0, Real tf, const State& x0, Real dt,
                              MatXd* J) {
    auto func = [=, this](const State& x) { return Propagate(odefunc, t0, tf, x, dt); };
    State xf = (J != nullptr) ? JacobianParallel(func, x0, *J) : func(x0);
    return xf;
  };

  IntegratorResult Integrator::PropagateEx(const ODE& odefunc, Real t0, Real tf, const VecX& x0,
                                           Real dt) {
    LUPNT_CHECK(dt > 0, "dt must be positive (magnitude)", "Integrator");  // treat as magnitude
    VecX x = x0;
    Real t = t0;
    int steps = 0;

    const double span = std::abs((tf - t0).val());
    if (span == 0) {
      // No propagation needed
      return {x0, t0, TerminationReason::ReachedTf, 0};
    }

    const double dir = (tf.val() > t0.val()) ? 1 : -1;  // forward or backward
    if (print_progress_) {
      ProgressBar pbar(100);
      pbar.SetDescription(dir > 0 ? "Integrating (forward)" : "Integrating (backward)");

      // initial early-termination check
      if (params_.terminate_if && params_.terminate_if(t, x)) {
        pbar.Finish();
        return {x, t, TerminationReason::UserCondition, steps};
      }

      while ((dir > 0 && t < tf) || (dir < 0 && t > tf)) {
        // signed step limited so we land exactly on tf
        Real h = dir * std::min(std::abs(dt.val()), std::abs((tf - t).val()));
        Real prev_h = h;

        x = Step(odefunc, t, State(x), h);
        t += prev_h;
        ++steps;

        if (params_.terminate_if && params_.terminate_if(t, x)) {
          // progress uses absolute distances in time
          int pct = int(std::round(std::abs((t - t0).val()) / span * 100.0));
          pbar.Update(std::min(100, std::max(0, pct)));
          pbar.Finish();
          return {x, t, TerminationReason::UserCondition, steps};
        }

        int pct = int(std::round(std::abs((t - t0).val()) / span * 100.0));
        pbar.Update(std::min(100, std::max(0, pct)));
      }
      pbar.Finish();
      return {x, t, TerminationReason::ReachedTf, steps};
    } else {
      // same logic without progress bar
      if (params_.terminate_if && params_.terminate_if(t, x)) {
        return {x, t, TerminationReason::UserCondition, steps};
      }

      while ((dir > 0 && t < tf) || (dir < 0 && t > tf)) {
        Real h = dir * std::min(std::abs(dt.val()), std::abs((tf - t).val()));
        Real prev_h = h;

        x = Step(odefunc, t, State(x), h);
        t += prev_h;
        ++steps;

        if (params_.terminate_if && params_.terminate_if(t, x)) {
          return {x, t, TerminationReason::UserCondition, steps};
        }
      }
      return {x, t, TerminationReason::ReachedTf, steps};
    }
  }

  IntegratorResult Integrator::PropagateEx(const ODE& f, Real t0, Real tf, const VecX& x0, Real dt,
                                           MatXd* J) {
    // One simple approach: propagate nominally and get Jacobian via your existing
    // path. If you need the Jacobian at the *actual stop time*, call PropagateEx
    // inside the Jacobian evaluator.
    auto func = [=, this](const VecX& x) { return PropagateEx(f, t0, tf, x, dt); };
    if (J) JacobianParallel([&](const VecX& x) { return func(x).x; }, x0, *J);
    return func(x0);
  }

  VecX Integrator::Propagate(const ODE& odefunc, Real t0, Real tf, const VecX& x0, Real dt) {
    IntegratorResult res = PropagateEx(odefunc, t0, tf, x0, dt);
    return res.x;
  };

  VecX Integrator::Propagate(const ODE& odefunc, Real t0, Real tf, const VecX& x0, Real dt,
                             MatXd* J) {
    auto func = [=, this](const VecX& x) { return Propagate(odefunc, t0, tf, x, dt); };
    VecX xf = (J != nullptr) ? JacobianParallel(func, x0, *J) : func(x0);
    return xf;
  };

  /**
   * @brief One step of Runge-Kutta Integration
   *
   * @param f  The ODE function to propagate
   * @param t  Time
   * @param x  The state to propagate
   * @param dt  Timestep
   */
  State RK4::Step(const ODE& f, Real t, const State& x, Real dt) {
    // usleep(1000);
    // return x + x;

    /* Evaluate `f` (i.e., `dx`) at the 4 locations defined bx the RK4 method */
    State k_1 = f(t, x) * dt;

    Real t1 = t + dt / 2.0;
    State x1 = x + k_1 / 2.0;
    State k_2 = f(t1, x1) * dt;

    Real t2 = t + dt / 2.0;
    State x2 = x + k_2 / 2.0;
    State k_3 = f(t2, x2) * dt;

    Real t3 = t + dt;
    State x3 = x + k_3;
    State k_4 = f(t3, x3) * dt;

    /* Average the 4 derivatives to approtimate `dx` */
    State dx = (k_1 + k_2 * 2.0 + k_3 * 2.0 + k_4) / 6.0;

    return x + dx;
  }

  /**
   * @brief One step of the 10-stage Shanks-type Runge-Kutta method.
   *
   * @note Despite the class name, this is a **7th order** method. Butcher's
   * barrier requires at least 11 stages for an explicit Runge-Kutta method to
   * attain order 8, and this one has 10, so order 7 is the best attainable at
   * this stage count -- the implementation is optimal, the name overstates it.
   * Measured in exact rational arithmetic the local truncation error is
   * O(h^8), i.e. global order 7; see cpp/test/numerics/test_integrator_order.cc.
   * The name is retained because IntegratorType::RK8 is public API and is
   * exposed in the Python bindings.
   *
   * @param f  The ODE function to propagate
   * @param t  Time
   * @param x  The state to propagate
   * @param dt  Timestep
   *
   * @ref
   * https://www.mathworks.com/matlabcentral/fileexchange/55431-runge-kutta-8th-order-integration
   */
  State RK8::Step(const ODE& f, Real t, const State& x, Real dt) {
    // 1
    State k_1 = f(t, x) * dt;

    // 2
    Real t1 = t + dt * (4.0 / 27.0);
    State x1 = x + k_1 * (4.0 / 27.0);
    State k_2 = f(t1, x1) * dt;

    // 3
    Real t2 = t + dt * (2.0 / 9);
    State x2 = x + (1.0 / 18) * (k_1 + 3 * k_2);
    State k_3 = f(t2, x2) * dt;

    // 4
    Real t3 = t + dt * (1.0 / 3);
    State x3 = x + (1.0 / 12) * (k_1 + 3 * k_3);
    State k_4 = f(t3, x3) * dt;

    // 5
    Real t4 = t + dt * (1.0 / 2);
    State x4 = x + (1.0 / 8) * (k_1 + 3 * k_4);
    State k_5 = f(t4, x4) * dt;

    // 6
    Real t5 = t + dt * (2.0 / 3);
    State x5 = x + (1.0 / 54) * (13 * k_1 - 27 * k_3 + 42 * k_4 + 8 * k_5);
    State k_6 = f(t5, x5) * dt;

    // 7
    Real t6 = t + dt * (1.0 / 6);
    State x6 = x + (1.0 / 4320) * (389 * k_1 - 54 * k_3 + 966 * k_4 - 824 * k_5 + 243 * k_6);
    State k_7 = f(t6, x6) * dt;

    // 8
    // The leading coefficient is -231, not -234: the row must sum to c_8 = 1
    // (the consistency condition c_i = sum_j a_ij), and
    // (-231 + 81 - 1164 + 656 - 122 + 800)/20 = 1 exactly, whereas -234 gives
    // 17/20 and silently drops the method to 3rd order. See
    // cpp/test/numerics/test_integrator_order.cc.
    Real t7 = t + dt;
    State x7
        = x + (1.0 / 20) * (-231 * k_1 + 81 * k_3 - 1164 * k_4 + 656 * k_5 - 122 * k_6 + 800 * k_7);
    State k_8 = f(t7, x7) * dt;

    // 9
    Real t8 = t + (5.0 / 6) * dt;
    State x8
        = x
          + (1.0 / 288)
                * (-127 * k_1 + 18 * k_3 - 678 * k_4 + 456 * k_5 - 9 * k_6 + 576 * k_7 + 4 * k_8);
    State k_9 = f(t8, x8) * dt;

    // 10
    Real t9 = t + dt;
    State x9 = x
               + (1.0 / 820)
                     * (1481 * k_1 - 81 * k_3 + 7104 * k_4 - 3376 * k_5 + 72 * k_6 - 5040 * k_7
                        - 60 * k_8 + 720 * k_9);
    State k_10 = f(t9, x9) * dt;

    /* Approtimate `dx` */
    State dx
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
  State IRKF::Step(const ODE& f, Real t, const State& x, Real dt) {
    int n = x.size();
    State x_new_low(n), x_new_high(n);

    int iter;
    for (iter = 0; iter < params_.max_iter; ++iter) {
      Update(f, t, x, dt, x_new_low, x_new_high);
      bool within_tolerance = ComputeRelError(x_new_low, x_new_high, dt);
      if (within_tolerance) break;
    }
    LUPNT_CHECK(iter < params_.max_iter,
                "IRKF did not converge, increase max_iter or reduce tolerance", "IRKF");

    return x_new_low;
  }

  bool IRKF::ComputeRelError(const State& x_new_low, const State& x_new_high, Real& dt) {
    // Compute the error norm
    double tol = 0.0;
    double error = 0.0;
    bool within_tolerance = true;
    double max_error = 0.0;

    // Non-conservative acceptance threshold. J.C. Butcher, Numerical Methods
    // for Ordinary Differential Equations, p291.
    //
    // The exponent must be evaluated in floating point: `order_` is an int, so
    // `(order_ + 1) / order_` is integer division and collapses to 1 for every
    // order >= 1, giving 5 instead of 5^1.25 ~ 7.48 for RKF45.
    double accept_thresh
        = std::pow(order_ + 1, static_cast<double>(order_ + 1) / static_cast<double>(order_));

    for (size_t i = 0; i < x_new_low.size(); ++i) {
      error = abs(x_new_high(i) - x_new_low(i));
      tol = max(params_.reltol * abs(x_new_high(i)), params_.abstol);
      max_error = std::max(max_error, error / tol);
      if ((error / tol) > accept_thresh) within_tolerance = false;
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

  void RKF45::Update(const ODE& f, Real t, const State& x, const Real dt, State& x_new_low,
                     State& x_new_high) {
    // Compute k1 to k6
    State k1 = f(t, x) * dt;
    State k2 = f(t + dt / 4.0, x + k1 * (1.0 / 4.0)) * dt;
    State k3 = f(t + dt * 3.0 / 8.0, x + k1 * (3.0 / 32.0) + k2 * (9.0 / 32.0)) * dt;
    State k4 = f(t + dt * 12.0 / 13.0,
                 x + k1 * (1932.0 / 2197.0) - k2 * (7200.0 / 2197.0) + k3 * (7296.0 / 2197.0))
               * dt;
    State k5 = f(t + dt, x + k1 * (439.0 / 216.0) - k2 * 8.0 + k3 * (3680.0 / 513.0)
                             - k4 * (845.0 / 4104.0))
               * dt;
    State k6 = f(t + dt / 2.0, x - k1 * (8.0 / 27.0) + k2 * 2.0 - k3 * (3544.0 / 2565.0)
                                   + k4 * (1859.0 / 4104.0) - k5 * (11.0 / 40.0))
               * dt;

    // 4th order solution
    x_new_low = x + k1 * (25.0 / 216.0) + k3 * (1408.0 / 2565.0) + k4 * (2197.0 / 4104.0)
                - k5 * (1.0 / 5.0);

    // 5th order solution
    x_new_high = x + k1 * (16.0 / 135.0) + k3 * (6656.0 / 12825.0) + k4 * (28561.0 / 56430.0)
                 - k5 * (9.0 / 50.0) + k6 * (2.0 / 55.0);
  }

  const std::array<std::array<double, 6>, 7> PD45::A_
      = {{{0, 0, 0, 0, 0, 0},
          {1.0 / 5, 0, 0, 0, 0, 0},
          {3.0 / 40, 9.0 / 40, 0, 0, 0, 0},
          {44.0 / 45, -56.0 / 15, 32.0 / 9, 0, 0, 0},
          {19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729, 0, 0},
          {9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656, 0},
          {35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84}}};

  const std::array<double, 7> PD45::b_
      = {35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84, 0};

  const std::array<double, 7> PD45::b_star_{
      5179.0 / 57600, 0, 7571.0 / 16695, 393.0 / 640, -92097.0 / 339200, 187.0 / 2100, 1.0 / 40};

  State PD45::Step(const ODE& f, Real t, const State& x, Real dt) {
    const int stages = 7;
    std::vector<State> k(stages);  // To store intermediate k values

    double safety = 0.9;   // Safety factor for step size adjustment
    double dt_min = 1e-3;  // Minimum allowable step size
    int iteration = 0;

    while (iteration < this->params_.max_iter) {
      // Compute intermediate stages
      k[0] = dt * f(t, x);
      for (int i = 1; i < stages; ++i) {
        State sum = State::Zero(x.size());
        // The stage time is the ROW SUM c_i = sum_j a_ij, not the first column
        // A_[i][0]. Deriving it here rather than storing a separate c_ array
        // keeps the two consistent by construction. Using A_[i][0] put stages
        // 4 and 5 at t + 2.95*dt and t + 2.85*dt -- outside the step -- which
        // is invisible for an autonomous RHS but drops the method to 1st order
        // as soon as f depends on t, as orbit dynamics with ephemeris terms
        // does. See cpp/test/numerics/test_integrator_order.cc.
        double c_i = 0.0;
        for (int j = 0; j < i; ++j) {
          sum += A_[i][j] * k[j];
          c_i += A_[i][j];
        }
        k[i] = dt * f(t + c_i * dt, x + sum);
      }

      // Compute the high-order and low-order solutions
      State y_high = x;  // High-order solution
      State y_low = x;   // Low-order solution
      for (int i = 0; i < stages; ++i) {
        y_high += b_[i] * k[i];
        y_low += b_star_[i] * k[i];
      }

      // Compute the error norm
      State error_vec = y_high - y_low;
      State scale
          = (Real(this->params_.abstol) + this->params_.reltol * x.cwiseAbs().array()).matrix();
      Real error_norm = (error_vec.cwiseQuotient(scale)).norm() / std::sqrt(x.size());

      // Adjust step size based on error norm
      if (error_norm <= 1.0) {
        // Step is accepted
        t += dt;
        return y_high;  // Return the high-order solution
      } else {
        // Step is rejected
        double factor = std::pow(safety / error_norm.val(), 0.25);
        dt *= std::max(factor, 0.1);  // Decrease step size
      }

      // Ensure the step size is not too small
      if (abs(dt.val()) < dt_min) {
        throw std::runtime_error("Step size too small, unable to proceed.");
      }

      iteration++;
    }

    throw std::runtime_error("Maximum iterations exceeded in PD45 step.");
  }

}  // namespace lupnt
