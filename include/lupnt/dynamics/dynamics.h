/**
 * @file dynamics.h
 * @author Stanford NAV LAB
 * @brief Interface for Dynamics
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include "lupnt/dynamics/propagator.h"
#include "lupnt/numerics/integrator.h"
#include "lupnt/physics/body.h"
#include "lupnt/physics/orbit_state.h"
#include "lupnt/physics/state.h"

namespace lupnt {

  // ****************************************************************************
  // Base Dyanmics Classes
  // ****************************************************************************

  // Dynamics Interface
  class IDynamics {
  public:
    virtual ~IDynamics() = default;

    // Interface
    virtual Ptr<IState> PropagateState(const Ptr<IState> &state, Real t0, Real tf,
                                       MatXd *stm = nullptr)
        = 0;
    virtual VecX Propagate(const VecX &x0, Real t0, Real tf, MatXd *stm = nullptr) = 0;
  };

  // Orbit Dynamics Interface
  class IOrbitDynamics : public IDynamics {
  public:
    virtual ~IOrbitDynamics() = default;

    // Overrides
    Ptr<IState> PropagateState(const Ptr<IState> &state, Real t0, Real tf,
                               MatXd *stm = nullptr) override;
    VecX Propagate(const VecX &x0, Real t0, Real tf, MatXd *stm = nullptr) override;

    // Interface
    virtual OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                                      Mat6d *stm = nullptr)
        = 0;
    virtual Vec6 Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm = nullptr) = 0;
    virtual MatX6 Propagate(const Vec6 &x0, Real t0, const VecX &tf, bool progress = false) = 0;

    // Implementations
    MatX6 Propagate(const MatX6 &x0, Real t0, Real tf);
  };

  // Analytical Orbit Dynamics Interface
  class IAnalyticalOrbitDynamics : public IOrbitDynamics {
  public:
    virtual ~IAnalyticalOrbitDynamics() = default;

    // Overrides
    using IOrbitDynamics::Propagate;
    MatX6 Propagate(const Vec6 &x0, Real t0, const VecX &tf, bool progress = false) override;

    // Interface
    virtual OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                                      Mat6d *stm = nullptr) override
        = 0;
    virtual Vec6 Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm = nullptr) override = 0;
  };

  class NumericalOrbitDynamics : public IOrbitDynamics {
  private:
    ODE odefunc_;
    NumericalPropagator propagator_;
    Real dt_ = 10.0;

  public:
    NumericalOrbitDynamics(ODE odefunc = nullptr, IntegratorType integ = default_integrator);
    void SetTimeStep(Real dt);
    Real GetTimeStep() const;
    void SetODEFunction(ODE odefunc);
    void SetIntegratorParams(IntegratorParams params) {
      propagator_.integrator->SetIntegratorParams(params);
    }

    // Overrides
    using IOrbitDynamics::Propagate;
    Vec6 Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm = nullptr) override;
    MatX6 Propagate(const Vec6 &x0, Real t0, const VecX &tf, bool progress = false) override;

    // Interface
    virtual Vec6 ComputeRates(Real t, const Vec6 &x) const = 0;
    virtual OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                                      Mat6d *stm = nullptr) override
        = 0;
  };

  // ****************************************************************************
  // Analytical Dynamics Classes
  // ****************************************************************************

  // Keplerian Dynamics
  class KeplerianDynamics : public IAnalyticalOrbitDynamics {
  private:
    Real GM_;

  public:
    KeplerianDynamics(Real GM);
    Vec6 PropagateClassicalOE(const Vec6 &coe, Real t0, Real tf, Mat6d *stm = nullptr);
    Vec6 PropagateQuasiNonsingOE(const Vec6 &qnsoe, Real t0, Real tf, Mat6d *stm = nullptr);
    Vec6 PropagateEquinoctialOE(const Vec6 &eqoe, Real t0, Real tf, Mat6d *stm = nullptr);

    using IAnalyticalOrbitDynamics::Propagate;
    Vec6 Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm = nullptr) override;
    OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                              Mat6d *stm = nullptr) override;
  };

  // Clohessy-Wiltshire Dynamics
  class ClohessyWiltshireDynamics : public IAnalyticalOrbitDynamics {
  private:
    Real a_, n_;
    Vec6 K_;
    Real t0_;

  public:
    ClohessyWiltshireDynamics(Real a, Real n);
    Mat6 ComputeMat(Real tf);

    using IAnalyticalOrbitDynamics::Propagate;
    Vec6 Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm = nullptr) override;
    OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                              Mat6d *stm = nullptr) override;
  };

  // Yamanaka-Ankersen Dynamics
  class YamanakaAnkersenDynamics : public IAnalyticalOrbitDynamics {
  private:
    Real a_, n_, e_, M0_;
    Vec6 K_;
    Real t0_;
    Vec6 rv_rtn_;

  public:
    YamanakaAnkersenDynamics(const ClassicalOE &coe_c, const CartesianOrbitState &rv_rtn, Real GM_);
    MatX ComputeMat(Real t);
    MatX ComputeInverseMat(Real t);

    using IAnalyticalOrbitDynamics::Propagate;
    Vec6 Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm = nullptr) override;
    OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                              Mat6d *stm = nullptr) override;
  };

  // Roe Geometric Mapping Dynamics
  class RoeGeometricMappingDynamics : public IAnalyticalOrbitDynamics {
  private:
    Real a_, e_, i_, w_, M0_;
    Real ex_, ey_, n_;
    Vec6 K_;
    Real t0_;

  public:
    RoeGeometricMappingDynamics(const ClassicalOE coe_c, const QuasiNonsingROE &roe, Real GM);
    MatX ComputeMat(Real t);

    using IAnalyticalOrbitDynamics::Propagate;
    Vec6 Propagate(const Vec6 &x0, Real t0, Real tf, Mat6d *stm = nullptr) override;
    OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                              Mat6d *stm = nullptr) override;
  };

  // ****************************************************************************
  // Numerical Dynamics Classes
  // ****************************************************************************

  // Cartesian Two-Body Dynamics
  class CartesianTwoBodyDynamics : public NumericalOrbitDynamics {
  private:
    Real GM_;

  public:
    CartesianTwoBodyDynamics(Real GM, IntegratorType integ = default_integrator);
    Vec6 ComputeRates(Real t, const Vec6 &x) const override;
    OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                              Mat6d *stm = nullptr) override;
  };

  // J2 Cartesian Two-Body Dynamics
  class J2CartTwoBodyDynamics : public NumericalOrbitDynamics {
  private:
    Real GM_, J2_, R_body_;

  public:
    J2CartTwoBodyDynamics(Real GM, Real J2, Real R_body, IntegratorType integ = default_integrator);
    Vec6 ComputeRates(Real t, const Vec6 &x) const override;
    OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                              Mat6d *stm = nullptr) override;
  };

  class J2KeplerianDynamics : public NumericalOrbitDynamics {
  private:
    Real GM_, J2_, R_body_;

  public:
    J2KeplerianDynamics(Real GM, Real J2, Real R_body, IntegratorType integ = default_integrator);
    Vec6 ComputeRates(Real t, const Vec6 &x) const override;
    OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                              Mat6d *stm = nullptr) override;
  };

  class MoonMeanDynamics : public NumericalOrbitDynamics {
  private:
    const double n3_ = 2.66e-6;
    const double J2_ = 2.03e-4;
    const double k_ = 0.98785;

  public:
    MoonMeanDynamics(IntegratorType integ = default_integrator);
    Vec6 ComputeRates(Real t, const Vec6 &x) const override;
    OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                              Mat6d *stm = nullptr) override;
  };

  struct NBodyDynamicsParams {
    Real mass;  // [kg] Spacecraft mass
    Real area;  // [m^2] Cross-sectional area
    Real CR;    // [-] Radiation pressure coefficient
    Real CD;    // [-] Drag coefficient
    bool use_srp = false;
    bool use_drag = false;
  };

  template <typename T = double> class NBodyDynamics : public NumericalOrbitDynamics {
  private:
    Frame frame_ = Frame::NONE;
    std::vector<BodyT<T>> bodies_;
    NumericalPropagator propagator;
    ODE odefunc;
    bool use_srp_ = false;
    bool use_drag_ = false;

    Real mass_;  // [kg] Spacecraft mass
    Real area_;  // [m^2] Cross-sectional area
    Real CR_;    // [-] Radiation pressure coefficient
    Real CD_;    // [-] Drag coefficient

  public:
    NBodyDynamics(IntegratorType integ = default_integrator);
    Vec6 ComputeRates(Real epoch, const Vec6 &x) const override;
    OrbitState PropagateState(const OrbitState &state, Real t0, Real tf,
                              Mat6d *stm = nullptr) override;

    void AddBody(const BodyT<T> &body) {
      for (auto &b : bodies_) {
        if (b.id == body.id) throw std::runtime_error("Body already added");
      }
      bodies_.push_back(body);
    }

    std::vector<BodyT<T>> GetBodies() { return bodies_; }

    void RemoveBody(const BodyT<T> &body) {
      for (auto it = bodies_.begin(); it != bodies_.end(); ++it) {
        if (it->id == body.id) {
          bodies_.erase(it);
          break;
        }
      }
    }

    void SetFrame(Frame frame) { frame_ = frame; }
    void GetFrame(Frame &frame) { frame = frame_; }

    void SetMass(Real mass) { mass_ = mass; }
    void SetArea(Real area) { area_ = area; }
    void SetSrpCoeff(Real CR) { CR_ = CR; }
    void SetDragCoeff(Real CD) { CD_ = CD; }

    void SetUseSrp(bool use_srp) { use_srp_ = use_srp; }
    void SetUseDrag(bool use_drag) { use_drag_ = use_drag; }
  };

}  // namespace lupnt
