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

  /**
   * @brief Interface for Dynamics
   *
   */
  class IDynamics {
  public:
    // without dt
    virtual ~IDynamics() = default;
    virtual void PropagateX(VecX &x, Real t0, Real tf) = 0;
    virtual void PropagateWithStmX(VecX &x, Real t0, Real tf, MatXd &stm) = 0;
  };

  /**
   * @brief Numerical Dynamics
   *
   */
  class NumericalDynamics : public IDynamics {
  protected:
    ODE odefunc_;
    NumericalPropagator propagator_;
    double dt_ = 0.0;

  public:
    NumericalDynamics(ODE odefunc, std::string integrator = "RK4")
        : odefunc_(odefunc), propagator_(integrator){};

    void SetTimeStep(Real dt) { dt_ = dt.val(); };
    void PropagateX(VecX &x, Real t0, Real tf);
    void PropagateWithStmX(VecX &x, Real t0, Real tf, MatXd &stm);

  protected:
    virtual VecX ComputeRates(Real t, const VecX &x) const = 0;
  };

  /********************************************
   * Analytical Orbit Dynamics
   **********************************************/

  /**
   * @brief Analytical Dynamics
   *
   */
  class AnalyticalDynamics : public IDynamics {
  public:
    virtual ~AnalyticalDynamics(){};
    virtual OrbitState CreateOrbitState(Vec6 &x) = 0;
    virtual void Propagate(OrbitState &state, Real t0, Real dt) = 0;
    virtual void PropagateWithSTM(OrbitState &state, Real t0, Real dt, Mat6d &stm) = 0;

    // Using fixed size vectors
    void Propagate(Vec6 &x, Real t0, Real dt) {
      OrbitState state = CreateOrbitState(x);
      Propagate(state, t0, dt);
      x = state.GetVec();
    }

    void PropagateWithSTM(Vec6 &x, Real t0, Real dt, Mat6d &stm) {
      OrbitState state = CreateOrbitState(x);
      PropagateWithSTM(state, t0, dt, stm);
      x = state.GetVec();
    }

    // arbitrary state size
    void PropagateX(VecX &x, Real t0, Real tf) {
      Vec6 x6 = x.head(6);
      Real dt = tf - t0;
      Propagate(x6, t0, dt);
      x.head(6) = x6;
    }

    void PropagateWithStmX(VecX &x, Real t0, Real tf, MatXd &stm) {
      Vec6 x6 = x.head(6);
      Real dt = tf - t0;
      Mat6d stm6;
      stm6 = stm.block(0, 0, 6, 6);
      PropagateWithSTM(x6, t0, dt, stm6);
      x.head(6) = x6;
      stm.block(0, 0, 6, 6) = stm6;
    }
  };

  /**
   * @brief Keplerian Dynamics (Todo: find a way to inherit this from
   * AnalyticalDynamics)
   *
   */
  class KeplerianDynamics {
  private:
    double mu_;

  public:
    KeplerianDynamics(double GM);

    // ClassicalOE
    void Propagate(ClassicalOE &state, Real dt);
    void PropagateWithStm(ClassicalOE &state, Real dt, Mat6d &stm);
    static Vec6 PropagateClassicalOE(Vec6 coe, Real dt, double GM);
    // QuasiNonsingOE
    void Propagate(QuasiNonsingOE &state, Real dt);
    void PropagateWithStm(QuasiNonsingOE &state, Real dt, Mat6d &stm);
    // EquinoctialOE
    void Propagate(EquinoctialOE &state, Real dt);
    void PropagateWithStm(EquinoctialOE &state, Real dt, Mat6d &stm);
  };

  class ClohessyWiltshireDynamics : public AnalyticalDynamics {
  private:
    Real a, n;
    VecX K;
    Real tInit;

  public:
    ClohessyWiltshireDynamics(Real a_in, Real n_in);
    void Propagate(OrbitState &state, Real tf);
    void Initialize(CartesianOrbitState &state, Real t0);
    MatX ComputeMat(Real t);
  };

  class YamanakaAnkersenDynamics : public AnalyticalDynamics {
  private:
    Real a, n, e, M0;
    VecX K;
    Real tInit;

  public:
    YamanakaAnkersenDynamics();
    void Propagate(CartesianOrbitState &state, Real tf);
    void Initialize(ClassicalOE &coe_c, CartesianOrbitState &rv_rtn, Real t0, double GM);
    MatX ComputeMat(Real t);
    MatX ComputeInverseMat(Real t);
  };

  class RoeGeometricMappingDynamics : public AnalyticalDynamics {
  private:
    Real a, e, i, w, M0, ex, ey, n;
    VecX K;
    Real tInit;

  public:
    RoeGeometricMappingDynamics();
    void Propagate(CartesianOrbitState &state, Real tf);
    void Initialize(ClassicalOE coe_c, QuasiNonsingROE &roe, Real t0, double GM);
    MatX ComputeMat(Real t);
  };

  /********************************************
   * Numerical Orbit Dynamics
   **********************************************/

  /**
   * @brief Numerical Dynamics
   *
   */
  class NumericalOrbitDynamics : public NumericalDynamics {
  private:
    OrbitStateRepres state_representation_;

  public:
    NumericalOrbitDynamics(ODE odefunc, OrbitStateRepres state_representation,
                           std::string integrator = "RK4")
        : NumericalDynamics(odefunc, integrator), state_representation_(state_representation) {}

    void Propagate(OrbitState &state, Real t0, Real tf, Real dt = 0.0);
    void Propagate(Vec6 &x, Real t0, Real tf, Real dt = 0.0);
    void PropagateWithStm(OrbitState &state, Real t0, Real tf, Real dt, Mat6d &stm);
    MatX Propagate(OrbitState &state, Real t0, VecX &tf, Real dt = 0.0, bool progress = false);
    MatX Propagate(Vec6 x, Real t0, VecX tf, Real dt = 0.0, bool progress = false);
    void PropagateWithStm(Vec6 &x, Real t0, Real tf, Real dt, Mat6d &stm);
    void PropagateWithStm(OrbitState &state, Real t0, Real tf, Mat6d &stm);
    void PropagateWithStm(Vec6 &x, Real t0, Real tf, Mat6d &stm);

    // with returns
    Vec6 PropagateR(Vec6 &x, Real t0, Real tf, Real dt);
    Vec6 PropagateWithStmR(Vec6 &x, Real t0, Real tf, Real dt, Mat6d &stm);
    Vec6 PropagateR(Vec6 &x, Real t0, Real tf);
    Vec6 PropagateWithStmR(Vec6 &x, Real t0, Real tf, Mat6d &stm);

  protected:
    virtual VecX ComputeRates(Real t, const VecX &x) const = 0;
  };

  class MoonFixedDynamics : public NumericalOrbitDynamics {
  private:
    double mu_;

  public:
    MoonFixedDynamics(double GM, std::string integrator = "RK4");
    VecX ComputeRates(Real t, const VecX &x) const;
  };

  class CartesianTwoBodyDynamics : public NumericalOrbitDynamics {
  private:
    double mu_;

  public:
    CartesianTwoBodyDynamics(double GM, std::string integrator = "RK4");
    VecX ComputeRates(Real t, const VecX &x) const;
  };

  class J2CartTwoBodyDynamics : public NumericalOrbitDynamics {
  private:
    double mu_, J2_, Rbody_;

  public:
    J2CartTwoBodyDynamics(double GM, double J2_in, double Rbody_in, std::string integrator = "RK4");
    VecX ComputeRates(Real t, const VecX &x) const;
  };

  class J2KeplerianDynamics : public NumericalOrbitDynamics {
  private:
    double GM, J2, Rbody;

  public:
    J2KeplerianDynamics(double GM, double J2_in, double Rbody_in, std::string integrator = "RK4");
    VecX ComputeRates(Real t, const VecX &x) const;
  };

  class MoonMeanDynamics : public NumericalOrbitDynamics {
  private:
    double n3 = 2.66e-6;
    double nM = 2.66e-6;
    double J2 = 2.03e-4;
    double k = 0.98785;

  public:
    MoonMeanDynamics(std::string integrator = "RK4");
    VecX ComputeRates(Real t, const VecX &x) const;
  };

  struct NBodyDynamicsParams {
    Real mass;  // [kg] Spacecraft mass
    Real area;  // [m^2] Cross-sectional area
    Real CR;    // [-] Radiation pressure coefficient
    Real CD;    // [-] Drag coefficient
    bool use_srp = false;
    bool use_drag = false;
  };

  class NBodyDynamics : public NumericalOrbitDynamics {
  private:
    Body central_body_;
    std::vector<Body> bodies_;
    NumericalPropagator propagator;
    ODE odefunc;
    bool use_srp_ = false;
    bool use_drag_ = false;

    Real mass_;  // [kg] Spacecraft mass
    Real area_;  // [m^2] Cross-sectional area
    Real CR_;    // [-] Radiation pressure coefficient
    Real CD_;    // [-] Drag coefficient

  public:
    NBodyDynamics(std::string integrator = "RK4");
    VecX ComputeRates(Real epoch, const VecX &x) const;

    void AddBody(const Body &body) {
      for (auto &b : bodies_) {
        assert(b.id != body.id && "Body already added");
      }
      bodies_.push_back(body);
    }

    void RemoveBody(const Body &body) {
      for (auto it = bodies_.begin(); it != bodies_.end(); ++it) {
        if (it->id == body.id) {
          bodies_.erase(it);
          break;
        }
      }
    }

    void SetPrimaryBody(const Body &body) {
      central_body_ = body;
      bodies_.push_back(body);
    }
    void SetMass(Real mass) { mass_ = mass; }
    void SetArea(Real area) { area_ = area; }
    void SetSrpCoeff(Real CR) { CR_ = CR; }
    void SetDragCoeff(Real CD) { CD_ = CD; }

    void SetUseSrp(bool use_srp) { use_srp_ = use_srp; }
    void SetUseDrag(bool use_drag) { use_drag_ = use_drag; }
  };

}  // namespace lupnt