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
  virtual void Propagate(VectorXreal &x, real t0, real tf) = 0;
  virtual void PropagateWithStm(VectorXreal &x, real t0, real tf,
                                MatrixXd &stm) = 0;
};

/**
 * @brief Interface for Orbit Dynamics (Todo: find a way to inherit this from
 * IDynamics)
 *
 */
class IOrbitDynamics {
 public:
  virtual ~IOrbitDynamics() = default;
  virtual void Propagate(OrbitState &state, real t0, real tf, real dt) = 0;
  virtual void Propagate(Vector6real &x, real t0, real tf, real dt) = 0;
  virtual void PropagateWithStm(OrbitState &state, real t0, real tf, real dt,
                                Matrix6d &stm) = 0;
  virtual void PropagateWithStm(Vector6real &x, real t0, real tf, real dt,
                                Matrix6d &stm) = 0;
  // without dt
  virtual void Propagate(OrbitState &state, real t0, real tf) = 0;
  virtual void Propagate(Vector6real &x, real t0, real tf) = 0;
  virtual void PropagateWithStm(OrbitState &state, real t0, real tf,
                                Matrix6d &stm) = 0;
  virtual void PropagateWithStm(Vector6real &x, real t0, real tf,
                                Matrix6d &stm) = 0;
};

class NumericalDynamics : public IOrbitDynamics {
 private:
  ODE odefunc_;
  double dt_;
  OrbitStateRepres state_representation_;
  NumericalPropagator propagator_;

 public:
  NumericalDynamics(ODE odefunc, OrbitStateRepres state_representation,
                    std::string integrator = "RK4")
      : state_representation_(state_representation),
        odefunc_(odefunc),
        propagator_(integrator){};

  // with dt
  void SetDt(double dt) { dt_ = dt; };
  void Propagate(OrbitState &state, real t0, real tf, real dt);
  void Propagate(Vector6real &x, real t0, real tf, real dt);
  void PropagateWithStm(OrbitState &state, real t0, real tf, real dt,
                        Matrix6d &stm);
  void PropagateWithStm(Vector6real &x, real t0, real tf, real dt,
                        Matrix6d &stm);

  // without dt (uses dt_)
  void Propagate(OrbitState &state, real t0, real tf);
  void Propagate(Vector6real &x, real t0, real tf);
  void PropagateWithStm(OrbitState &state, real t0, real tf, Matrix6d &stm);
  void PropagateWithStm(Vector6real &x, real t0, real tf, Matrix6d &stm);

 protected:
  virtual VectorXreal ComputeRates(real t, const VectorXreal &x) const = 0;
};

class IAnalyticalDynamics {
 public:
  virtual ~IAnalyticalDynamics(){};
  virtual void Propagate(OrbitState &state, real t0, real dt) = 0;
  virtual void Propagate(Vector6real &x, real t0, real dt) = 0;
  virtual void PropagateWithSTM(OrbitState &state, real t0, real dt) = 0;
  virtual void PropagateWithSTM(Vector6real &x, real t0, real dt) = 0;
};

/**
 * @brief Keplerian Dynamics (Todo: find a way to inherit this from
 * IAnalyticalDynamics)
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

class MoonFixedDynamics : public NumericalDynamics {
 private:
  double mu_;

 public:
  MoonFixedDynamics(double mu, std::string integrator = "RK4");
  VectorXreal ComputeRates(real t, const VectorXreal &x) const;
};

class CartesianTwoBodyDynamics : public NumericalDynamics {
 private:
  double mu_;

 public:
  CartesianTwoBodyDynamics(double mu, std::string integrator = "RK4");
  VectorXreal ComputeRates(real t, const VectorXreal &x) const;
};

class J2CartesianTwoBodyDynamics : public NumericalDynamics {
 private:
  double mu_, J2_, Rbody_;

 public:
  J2CartesianTwoBodyDynamics(double mu, double J2_in, double Rbody_in,
                             std::string integrator = "RK4");
  VectorXreal ComputeRates(real t, const VectorXreal &x) const;
};

class J2KeplerianDynamics : public NumericalDynamics {
 private:
  double mu, J2, Rbody;

 public:
  J2KeplerianDynamics(double mu, double J2_in, double Rbody_in,
                      std::string integrator = "RK4");
  VectorXreal ComputeRates(real t, const VectorXreal &x) const;
};

class MoonMeanDynamics : public NumericalDynamics {
 private:
  double n3 = 2.66e-6;
  double nM = 2.66e-6;
  double J2 = 2.03e-4;
  double k = 0.98785;

 public:
  MoonMeanDynamics(std::string integrator = "RK4");
  VectorXreal ComputeRates(real t, const VectorXreal &x) const;
};

class ClohessyWiltshireDynamics : public IAnalyticalDynamics {
 private:
  real a, n;
  VectorXreal K;
  real tInit;

 public:
  ClohessyWiltshireDynamics(real a_in, real n_in);
  void Propagate(OrbitState &state, real tf);
  void Initialize(CartesianOrbitState &state, real t0);
  MatrixXreal ComputeMatrix(real t);
};

class YamanakaAnkersenDynamics : public IAnalyticalDynamics {
 private:
  real a, n, e, M0;
  VectorXreal K;
  real tInit;

 public:
  YamanakaAnkersenDynamics();
  void Propagate(CartesianOrbitState &state, real tf);
  void Initialize(ClassicalOE &coe_c, CartesianOrbitState &rv_rtn, real t0,
                  double mu);
  MatrixXreal ComputeMatrix(real t);
  MatrixXreal ComputeInverseMatrix(real t);
};

class RoeGeometricMappingDynamics : public IAnalyticalDynamics {
 private:
  real a, e, i, w, M0, ex, ey, n;
  VectorXreal K;
  real tInit;

 public:
  RoeGeometricMappingDynamics();
  void Propagate(CartesianOrbitState &state, real tf);
  void Initialize(ClassicalOE coe_c, QuasiNonsingularROE &roe, real t0,
                  double mu);
  MatrixXreal ComputeMatrix(real t);
};

struct Body {
  std::string name;
  double mu;
  double R;
  BodyId id;
  bool sphericalHarmonics;
  bool normalized;
  int n_max;
  int m_max;
  MatrixXd Cnm;
  MatrixXd Snm;

  static Body Moon(int n_max = 0, int m_max = 0);
  static Body Earth(int n_max = 0, int m_max = 0);
};

class NBodyDynamics : public NumericalDynamics {
 private:
  Body centralBody;
  std::vector<Body> bodies;
  NumericalPropagator propagator;
  ODE odefunc;

 public:
  NBodyDynamics(std::string integrator = "RK4");
  VectorXreal ComputeRates(real t, const VectorXreal &x) const;

  void AddBody(const Body &body);
  void SetCentralBody(const Body &body);

  Vector3real ComputeNBodyGravity(const real t, const VectorXreal &rv) const;
  Vector3real ComputeSolarRadiationPressure(const Vector3real &r_body2sc,
                                            const Vector3real &r_sun2sc,
                                            double R_body, double R_SUN,
                                            double m, double CR,
                                            double area) const;
};

}  // namespace lupnt
