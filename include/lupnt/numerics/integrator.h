/**
 * @file integrator.h
 * @author Stanford NAV LAB
 * @brief Integrator interfaces
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <functional>

#include "lupnt/core/constants.h"

namespace lupnt {

  using ODE = std::function<VecX(Real, const VecX&)>;

  enum class IntegratorType {
    RK4,
    RK8,
    RKF45,
  };
  constexpr IntegratorType default_integrator = IntegratorType::RK4;

  class IntegratorParams {
  public:
    int max_iter = 20;
    double abstol = 1e-6;
    double reltol = 1e-6;

    IntegratorParams() = default;
    IntegratorParams(int max_iter, double abstol, double reltol)
        : max_iter(max_iter), abstol(abstol), reltol(reltol) {
      CheckIntegratorParams();
    };
    void CheckIntegratorParams();
  };

  class IIntegrator {
  protected:
    IntegratorParams params_;

  public:
    virtual VecX Step(const ODE& f, Real t, const VecX& x, Real dt) = 0;
    void SetIntegratorParams(IntegratorParams params) { params_ = params; };
    virtual ~IIntegrator() {};
  };

  // Runge-Kutta Integrators
  class RK4 : public IIntegrator {
  public:
    VecX Step(const ODE& f, Real t, const VecX& x, Real dt);
  };

  class RK8 : public IIntegrator {
  public:
    VecX Step(const ODE& f, Real t, const VecX& x, Real dt);
  };

  // Runge-Kutta-Fehlberg Integrators with adaptive step size
  class IRKF : public IIntegrator {
  private:
    int order_;

  public:
    IRKF() = default;
    IRKF(IntegratorParams params, int order) : order_(order) { SetIntegratorParams(params); };
    VecX Step(const ODE& f, Real t, const VecX& x, Real dt) override;
    bool ComputeRelError(const VecX& x_new_low, const VecX& x_new_high, Real dt);
    virtual void Update(const ODE& f, Real t, const VecX& x, Real dt, VecX& x_new_low,
                        VecX& x_new_high)
        = 0;
    virtual ~IRKF() = default;
  };

  class RKF45 : public IRKF {
  public:
    RKF45(IntegratorParams params) : IRKF(params, 4) {};
    void Update(const ODE& f, Real t, const VecX& x, Real dt, VecX& x_new_low,
                VecX& x_new_high) override;
  };

}  // namespace lupnt
