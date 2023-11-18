/**
 * @file Dynamics.h
 * @author Stanford NAV LAB
 * @brief Interface for Dynamics
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <autodiff/forward/real.hpp>

#include "lupnt/dynamics/Propagator.h"
#include "lupnt/dynamics/gravity_field.h"
#include "lupnt/numerics/integrator.h"
#include "lupnt/physics/orbit_state.h"
#include "lupnt/physics/state.h"

namespace ad = autodiff;

namespace lupnt {

/**
 * @brief Interface for Dynamics
 *
 */
class IDynamics {
 public:
  // without dt
  virtual void Propagate(ad::VectorXreal &x, ad::real t0, ad::real tf) = 0;
  virtual void PropagateWithStm(ad::VectorXreal &x, ad::real t0, ad::real tf,
                                Eigen::MatrixXd &stm) = 0;
};

/**
 * @brief Interface for Orbit Dynamics (Todo: find a way to inherit this from
 * IDynamics)
 *
 */
class IOrbitDynamics {
 public:
  virtual ~IOrbitDynamics() = default;
  virtual void Propagate(OrbitState &state, ad::real t0, ad::real tf,
                         ad::real dt) = 0;
  virtual void Propagate(ad::Vector6real &x, ad::real t0, ad::real tf,
                         ad::real dt) = 0;
  virtual void PropagateWithStm(OrbitState &state, ad::real t0, ad::real tf,
                                ad::real dt, Eigen::Matrix6d &stm) = 0;
  virtual void PropagateWithStm(ad::Vector6real &x, ad::real t0, ad::real tf,
                                ad::real dt, Eigen::Matrix6d &stm) = 0;
  // without dt
  virtual void Propagate(OrbitState &state, ad::real t0, ad::real tf) = 0;
  virtual void Propagate(ad::Vector6real &x, ad::real t0, ad::real tf) = 0;
  virtual void PropagateWithStm(OrbitState &state, ad::real t0, ad::real tf,
                                Eigen::Matrix6d &stm) = 0;
  virtual void PropagateWithStm(ad::Vector6real &x, ad::real t0, ad::real tf,
                                Eigen::Matrix6d &stm) = 0;
};

class NumericalDynamics : public IOrbitDynamics {
 private:
  ODE odefunc_;
  double dt_;
  OrbitStateRepres state_representation_;
  NumericalPropagator integrator_;

 public:
  NumericalDynamics(ODE odefunc, OrbitStateRepres state_representation,
                    std::string integrator = "RK4")
      : state_representation_(state_representation),
        odefunc_(odefunc),
        integrator_(integrator){};

  // with dt
  void SetDt(double dt) { dt_ = dt; };
  void Propagate(OrbitState &state, ad::real t0, ad::real tf, ad::real dt);
  void Propagate(ad::Vector6real &x, ad::real t0, ad::real tf, ad::real dt);
  void PropagateWithStm(OrbitState &state, ad::real t0, ad::real tf,
                        ad::real dt, Eigen::Matrix6d &stm);
  void PropagateWithStm(ad::Vector6real &x, ad::real t0, ad::real tf,
                        ad::real dt, Eigen::Matrix6d &stm);

  // without dt (uses dt_)
  void Propagate(OrbitState &state, ad::real t0, ad::real tf);
  void Propagate(ad::Vector6real &x, ad::real t0, ad::real tf);
  void PropagateWithStm(OrbitState &state, ad::real t0, ad::real tf,
                        Eigen::Matrix6d &stm);
  void PropagateWithStm(ad::Vector6real &x, ad::real t0, ad::real tf,
                        Eigen::Matrix6d &stm);

 protected:
  virtual ad::VectorXreal ComputeRates(ad::real t,
                                       const ad::VectorXreal &x) const = 0;
};

class IAnalyticalDynamics {
 public:
  virtual ~IAnalyticalDynamics(){};
  virtual void Propagate(OrbitState &state, ad::real t0, ad::real dt) = 0;
  virtual void Propagate(ad::Vector6real &x, ad::real t0, ad::real dt) = 0;
  virtual void PropagateWithSTM(OrbitState &state, ad::real t0,
                                ad::real dt) = 0;
  virtual void PropagateWithSTM(ad::Vector6real &x, ad::real t0,
                                ad::real dt) = 0;
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
  void Propagate(ClassicalOE &state, ad::real dt);
  void PropagateWithStm(ClassicalOE &state, ad::real dt, Eigen::Matrix6d &stm);
  // QuasiNonsingularOE
  void Propagate(QuasiNonsingularOE &state, ad::real dt);
  void PropagateWithStm(QuasiNonsingularOE &state, ad::real dt,
                        Eigen::Matrix6d &stm);
  // NonsingularOE
  void Propagate(NonsingularOE &state, ad::real dt);
  void PropagateWithStm(NonsingularOE &state, ad::real dt,
                        Eigen::Matrix6d &stm);
  // EquinoctialOE
  void Propagate(EquinoctialOE &state, ad::real dt);
  void PropagateWithStm(EquinoctialOE &state, ad::real dt,
                        Eigen::Matrix6d &stm);
};

class MoonFixedDynamics : public NumericalDynamics {
 private:
  double mu_;

 public:
  MoonFixedDynamics(double mu, std::string integrator = "RK4");
  ad::VectorXreal ComputeRates(ad::real t, const ad::VectorXreal &x) const;
};

class CartesianTwoBodyDynamics : public NumericalDynamics {
 private:
  double mu_;

 public:
  CartesianTwoBodyDynamics(double mu, std::string integrator = "RK4");
  ad::VectorXreal ComputeRates(ad::real t, const ad::VectorXreal &x) const;
};

class J2CartesianTwoBodyDynamics : public NumericalDynamics {
 private:
  double mu, J2, Rbody;

 public:
  J2CartesianTwoBodyDynamics(double mu, double J2_in, double Rbody_in,
                             std::string integrator = "RK4");
  ad::VectorXreal ComputeRates(ad::real t, const ad::VectorXreal &x) const;
};

class J2KeplerianDynamics : public NumericalDynamics {
 private:
  double mu, J2, Rbody;

 public:
  J2KeplerianDynamics(double mu, double J2_in, double Rbody_in,
                      std::string integrator = "RK4");
  ad::VectorXreal ComputeRates(ad::real t, const ad::VectorXreal &x) const;
};

class MoonMeanDynamics : public NumericalDynamics {
 private:
  double n3 = 2.66e-6;
  double nM = 2.66e-6;
  double J2 = 2.03e-4;
  double k = 0.98785;

 public:
  MoonMeanDynamics(std::string integrator = "RK4");
  ad::VectorXreal ComputeRates(ad::real t, const ad::VectorXreal &x) const;
};

class ClohessyWiltshireDynamics : public IAnalyticalDynamics {
 private:
  ad::real a, n;
  ad::VectorXreal K;
  ad::real tInit;

 public:
  ClohessyWiltshireDynamics(ad::real a_in, ad::real n_in);
  void Propagate(OrbitState &state, ad::real tf);
  void Initialize(CartesianOrbitState &state, ad::real t0);
  ad::MatrixXreal ComputeMatrix(ad::real t);
};

class YamanakaAnkersenDynamics : public IAnalyticalDynamics {
 private:
  ad::real a, n, e, M0;
  ad::VectorXreal K;
  ad::real tInit;

 public:
  YamanakaAnkersenDynamics();
  void Propagate(CartesianOrbitState &state, ad::real tf);
  void Initialize(ClassicalOE &coe_c, CartesianOrbitState &rv_rtn, ad::real t0,
                  double mu);
  ad::MatrixXreal ComputeMatrix(ad::real t);
  ad::MatrixXreal ComputeInverseMatrix(ad::real t);
};

class RoeGeometricMappingDynamics : public IAnalyticalDynamics {
 private:
  ad::real a, e, i, w, M0, ex, ey, n;
  ad::VectorXreal K;
  ad::real tInit;

 public:
  RoeGeometricMappingDynamics();
  void Propagate(CartesianOrbitState &state, ad::real tf);
  void Initialize(ClassicalOE coe_c, QuasiNonsingularROE &roe, ad::real t0,
                  double mu);
  ad::MatrixXreal ComputeMatrix(ad::real t);
};

struct Body {
  std::string name;
  double mu_;
  double R;
  BodyId id;
  bool sphericalHarmonics;
  bool normalized;
  int n_max;
  int m_max;
  Eigen::MatrixXd Cnm;
  Eigen::MatrixXd Snm;

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
  ad::VectorXreal ComputeRates(ad::real t, const ad::VectorXreal &x) const;

  void AddBody(const Body &body);
  void SetCentralBody(const Body &body);

  ad::Vector3real ComputeNBodyGravity(const ad::real t,
                                      const ad::VectorXreal &rv) const;
  ad::Vector3real ComputeSolarRadiationPressure(
      const ad::Vector3real &r_body2sc, const ad::Vector3real &r_sun2sc,
      double R_body, double R_SUN, double m, double CR, double area) const;
};

}  // namespace lupnt
