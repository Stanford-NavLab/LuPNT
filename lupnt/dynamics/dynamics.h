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

#include "lupnt/dynamics/gravity_field.h"
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
  virtual void PropagateX(VectorX &x, real t0, real tf) = 0;
  virtual void PropagateWithStmX(VectorX &x, real t0, real tf,
                                 MatrixXd &stm) = 0;
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

  void SetTimeStep(real dt) { dt_ = dt.val(); };
  void PropagateX(VectorX &x, real t0, real tf);
  void PropagateWithStmX(VectorX &x, real t0, real tf, MatrixXd &stm);

 protected:
  virtual VectorX ComputeRates(real t, const VectorX &x) const = 0;
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
  virtual OrbitState CreateOrbitState(Vector6 &x) = 0;
  virtual void Propagate(OrbitState &state, real t0, real dt) = 0;
  virtual void PropagateWithSTM(OrbitState &state, real t0, real dt,
                                Matrix6d &stm) = 0;

  // Using fixed size vectors
  void Propagate(Vector6 &x, real t0, real dt) {
    OrbitState state = CreateOrbitState(x);
    Propagate(state, t0, dt);
    x = state.GetVector();
  }

  void PropagateWithSTM(Vector6 &x, real t0, real dt, Matrix6d &stm) {
    OrbitState state = CreateOrbitState(x);
    PropagateWithSTM(state, t0, dt, stm);
    x = state.GetVector();
  }

  // arbitrary state size
  void PropagateX(VectorX &x, real t0, real tf) {
    Vector6 x6 = x.head(6);
    real dt = tf - t0;
    Propagate(x6, t0, dt);
    x.head(6) = x6;
  }

  void PropagateWithStmX(VectorX &x, real t0, real tf, MatrixXd &stm) {
    Vector6 x6 = x.head(6);
    real dt = tf - t0;
    Matrix6d stm6;
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
  KeplerianDynamics(double mu);

  // ClassicalOE
  void Propagate(ClassicalOE &state, real dt);
  void PropagateWithStm(ClassicalOE &state, real dt, Matrix6d &stm);
  // QuasiNonsingularOE
  void Propagate(QuasiNonsingularOE &state, real dt);
  void PropagateWithStm(QuasiNonsingularOE &state, real dt, Matrix6d &stm);
  // EquinoctialOE
  void Propagate(EquinoctialOE &state, real dt);
  void PropagateWithStm(EquinoctialOE &state, real dt, Matrix6d &stm);
};

class ClohessyWiltshireDynamics : public AnalyticalDynamics {
 private:
  real a, n;
  VectorX K;
  real tInit;

 public:
  ClohessyWiltshireDynamics(real a_in, real n_in);
  void Propagate(OrbitState &state, real tf);
  void Initialize(CartesianOrbitState &state, real t0);
  MatrixX ComputeMatrix(real t);
};

class YamanakaAnkersenDynamics : public AnalyticalDynamics {
 private:
  real a, n, e, M0;
  VectorX K;
  real tInit;

 public:
  YamanakaAnkersenDynamics();
  void Propagate(CartesianOrbitState &state, real tf);
  void Initialize(ClassicalOE &coe_c, CartesianOrbitState &rv_rtn, real t0,
                  double mu);
  MatrixX ComputeMatrix(real t);
  MatrixX ComputeInverseMatrix(real t);
};

class RoeGeometricMappingDynamics : public AnalyticalDynamics {
 private:
  real a, e, i, w, M0, ex, ey, n;
  VectorX K;
  real tInit;

 public:
  RoeGeometricMappingDynamics();
  void Propagate(CartesianOrbitState &state, real tf);
  void Initialize(ClassicalOE coe_c, QuasiNonsingularROE &roe, real t0,
                  double mu);
  MatrixX ComputeMatrix(real t);
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
      : NumericalDynamics(odefunc, integrator),
        state_representation_(state_representation) {}

  // with dt
  void Propagate(OrbitState &state, real t0, real tf, real dt);
  void Propagate(Vector6 &x, real t0, real tf, real dt);
  void PropagateWithStm(OrbitState &state, real t0, real tf, real dt,
                        Matrix6d &stm);
  void PropagateWithStm(Vector6 &x, real t0, real tf, real dt, Matrix6d &stm);

  // without dt (uses dt_)
  void Propagate(OrbitState &state, real t0, real tf);
  void Propagate(Vector6 &x, real t0, real tf);
  MatrixX Propagate(OrbitState &state, real t0, VectorX &tf);
  MatrixX Propagate(Vector6 &x, real t0, VectorX &tf);
  void PropagateWithStm(OrbitState &state, real t0, real tf, Matrix6d &stm);
  void PropagateWithStm(Vector6 &x, real t0, real tf, Matrix6d &stm);

 protected:
  virtual VectorX ComputeRates(real t, const VectorX &x) const = 0;
};

class MoonFixedDynamics : public NumericalOrbitDynamics {
 private:
  double mu_;

 public:
  MoonFixedDynamics(double mu, std::string integrator = "RK4");
  VectorX ComputeRates(real t, const VectorX &x) const;
};

class CartesianTwoBodyDynamics : public NumericalOrbitDynamics {
 private:
  double mu_;

 public:
  CartesianTwoBodyDynamics(double mu, std::string integrator = "RK4");
  VectorX ComputeRates(real t, const VectorX &x) const;
};

class J2CartesianTwoBodyDynamics : public NumericalOrbitDynamics {
 private:
  double mu_, J2_, Rbody_;

 public:
  J2CartesianTwoBodyDynamics(double mu, double J2_in, double Rbody_in,
                             std::string integrator = "RK4");
  VectorX ComputeRates(real t, const VectorX &x) const;
};

class J2KeplerianDynamics : public NumericalOrbitDynamics {
 private:
  double mu, J2, Rbody;

 public:
  J2KeplerianDynamics(double mu, double J2_in, double Rbody_in,
                      std::string integrator = "RK4");
  VectorX ComputeRates(real t, const VectorX &x) const;
};

class MoonMeanDynamics : public NumericalOrbitDynamics {
 private:
  double n3 = 2.66e-6;
  double nM = 2.66e-6;
  double J2 = 2.03e-4;
  double k = 0.98785;

 public:
  MoonMeanDynamics(std::string integrator = "RK4");
  VectorX ComputeRates(real t, const VectorX &x) const;
};

class NBodyDynamics : public NumericalOrbitDynamics {
 private:
  Body central_body_;
  std::vector<Body> bodies_;
  NumericalPropagator propagator;
  ODE odefunc;
  bool use_srp_ = false;
  double mass_ = 100.0;      // s/c mass [kg]
  double CR_ = 1.5;          // solar radiation pressure coefficient
  double area_ = 1.0 / 1e6;  // solar radiation pressure area [km^2]

 public:
  NBodyDynamics(std::string integrator = "RK4");
  VectorX ComputeRates(real epoch, const VectorX &x) const;

  void AddBody(const Body &body) {
    for (auto &b : bodies_) {
      assert(b.id != body.id && "Body already added");
    }
    bodies_.push_back(body);
  }
  void SetPrimaryBody(const Body &body) {
    central_body_ = body;
    bodies_.push_back(body);
  }
  void SetMass(double mass) { mass_ = mass; }
  void SetSrpArea(double area) { area_ = area; }
  void SetSrpCoeff(double CR) { CR_ = CR; }
  void SetSolarRadiationPressure(bool use_srp) { use_srp_ = use_srp; }

  Vector3 ComputeNBodyGravity(const real epoch, const VectorX &rv) const;
  Vector3 ComputeSolarRadiationPressure(const Vector3 &r_body2sc,
                                        const Vector3 &r_sun2sc,
                                        const real B_srp, double R_body) const;
};

}  // namespace lupnt
