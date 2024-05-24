#include <lupnt/core/constants.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/physics/orbit_state.h>
#include <lupnt/physics/orbit_state_utils.h>

#include <catch2/catch_test_macros.hpp>
#include <iostream>

#include "test_utils.cc"

using namespace lupnt;

// Keplerian Dynamics with Classical Orbital Elements
TEST_CASE("Test_KeplerianDynamics_ClassicalOE") {
  // Classical orbital elements
  real a = 6541.4;                  // [km]
  real e = 0.6;                     // [-]
  real i = 65.5 * RAD_PER_DEG;      // [rad]
  real Omega = 90.0 * RAD_PER_DEG;  // [rad]
  real w = 0.0 * RAD_PER_DEG;       // [rad]
  real M = 0.0 * RAD_PER_DEG;       // [rad]

  double mu = MU_MOON;
  ClassicalOE coe_state({a, e, i, Omega, w, M}, Frame::MI);
  Vector6 coe_analytical = coe_state.GetVector();

  // Keplerian dynamics
  auto kep_dyn = KeplerianDynamics(mu);

  // Propagation
  real dt = 10.0;                   // [s]
  real n = sqrt(mu / pow(a, 3.0));  // [rad/s]
  for (int i = 0; i < 100; i++) {
    kep_dyn.Propagate(coe_state, dt);
    coe_analytical(5) = wrapToPi(coe_analytical(5) + n * dt);

    EXPECT_NEAR_ADVEC(coe_state.GetVector(), coe_analytical, 1e-6);
    dt += 2.0;
  }

  // Propagation with STM
  auto propagate_function = [&](VectorX &vec, real dt) {
    ClassicalOE state(vec, Frame::MI);
    kep_dyn.Propagate(state, dt);
    vec = state.GetVector();
  };

  // Propagation with STM
  Matrix6d stm;
  Matrix6d stm_numerical;
  for (int i = 0; i < 100; i++) {
    NumericalJacobian(propagate_function, coe_state.GetVector(), dt,
                      stm_numerical, 1e-6);
    kep_dyn.PropagateWithStm(coe_state, dt, stm);

    EXPECT_NEAR_EIGENMAT(stm, stm_numerical, 1e-6);
    dt += 2.0;
  }
}

// CartesianTwoBodyDynamics
TEST_CASE("Test_CartesianTwoBodyDynamics") {
  // Classical orbital elements
  real a = 6541.4;                  // [km]
  real e = 0.6;                     // [-]
  real i = 65.5 * RAD_PER_DEG;      // [rad]
  real Omega = 90.0 * RAD_PER_DEG;  // [rad]
  real w = 0.0 * RAD_PER_DEG;       // [rad]
  real M = 0.0 * RAD_PER_DEG;       // [rad]

  double mu = MU_MOON;
  ClassicalOE coe_state({a, e, i, Omega, w, M}, Frame::MI);
  CartesianOrbitState cart_state = ClassicalToCartesian(coe_state, mu);
  Vector6 cart_vector = cart_state.GetVector();
  VectorX cart_vector_kep;

  // Two body dynamics
  auto kep_dyn = KeplerianDynamics(mu);
  auto tb_dyn = CartesianTwoBodyDynamics(mu);

  // Propagation
  real dt = 10.0;  // [s]
  for (int i = 0; i < 5; i++) {
    kep_dyn.Propagate(coe_state, dt);
    tb_dyn.Propagate(cart_state, 0.0, dt, 1.0);
    tb_dyn.Propagate(cart_vector, 0.0, dt, 1.0);

    cart_vector_kep = ClassicalToCartesian(coe_state.GetVector(), mu);
    EXPECT_NEAR_ADVEC(cart_vector_kep, cart_state.GetVector(), 1e-6);
    EXPECT_NEAR_ADVEC(cart_vector_kep, cart_vector, 1e-6);
  }

  // Propagation with STM
  auto propagate_function = [&](VectorX &vec, real dt) {
    CartesianOrbitState state(vec, Frame::MI);
    tb_dyn.Propagate(state, 0.0, dt, 0.1);
    vec = state.GetVector();
  };

  // Propagation with STM
  Matrix6d stm_state, stm_vector, stm_numerical;
  Vector6 cart_numerical = cart_vector;

  for (int i = 0; i < 5; i++) {
    NumericalJacobian(propagate_function, cart_numerical, dt, stm_numerical);

    tb_dyn.Propagate(cart_numerical, 0.0, dt, 0.1);
    tb_dyn.PropagateWithStm(cart_state, 0.0, dt, 0.1, stm_state);
    tb_dyn.PropagateWithStm(cart_vector, 0.0, dt, 0.1, stm_vector);

    EXPECT_NEAR_ADVEC(cart_numerical, cart_state.GetVector(), 1e-6);
    EXPECT_NEAR_ADVEC(cart_numerical, cart_vector, 1e-6);

    EXPECT_NEAR_EIGENMAT(stm_numerical, stm_state, 1e-5);
    EXPECT_NEAR_EIGENMAT(stm_numerical, stm_vector, 1e-5);
  }
}

/**
// J2 Cartesian dynamics
TEST_CASE("Test_CartesianJ2Dynamics") {
  // Classical orbital elements
  real a = 6541.4;                  // [km]
  real e = 0.6;                     // [-]
  real i = 65.5 * RAD_PER_DEG;      // [rad]
  real Omega = 90.0 * RAD_PER_DEG;  // [rad]
  real w = 0.0 * RAD_PER_DEG;       // [rad]
  real M = 0.0 * RAD_PER_DEG;       // [rad]

  double mu = MU_MOON;
  double J2 = J2_MOON;
  double Rbody = R_MOON;
  real dt = 10.0;

  ClassicalOE coe_state({a, e, i, Omega, w, M}, CoordSystem::MI);
  CartesianOrbitState cart_state = ClassicalToCartesian(coe_state, mu);
  Vector6 cart_vector = cart_state.GetVector();
  VectorX cart_vector_kep;

  // Two body dynamics
  auto j2_dyn = J2CartesianTwoBodyDynamics(mu, J2, Rbody, "RK4");
  auto j2_kep_dyn = J2KeplerianDynamics(mu, J2, Rbody, "RK4");

  // comparison
  for (int i = 0; i < 5; i++) {
    j2_kep_dyn.Propagate(coe_state, 0.0, dt, 1.0);
    j2_dyn.Propagate(cart_state, 0.0, dt, 1.0);
    j2_dyn.Propagate(cart_vector, 0.0, dt, 1.0);

    cart_vector_kep = ClassicalToCartesian(coe_state.GetVector(), mu);
    EXPECT_NEAR_ADVEC(cart_vector_kep, cart_state.GetVector(), 1e-6);
    EXPECT_NEAR_ADVEC(cart_vector_kep, cart_vector, 1e-6);
  }

  // Propagation with STM
  auto propagate_function = [&](VectorX &vec, real dt) {
    CartesianOrbitState state(vec, CoordSystem::MI);
    j2_dyn.Propagate(state, 0.0, dt, 0.1);
    vec = state.GetVector();
  };

  // Propagation with STM
  Matrix6d stm_state, stm_vector, stm_numerical;
  Vector6 cart_numerical = cart_vector;

  for (int i = 0; i < 5; i++) {
    NumericalJacobian(propagate_function, cart_numerical, dt, stm_numerical);

    j2_dyn.Propagate(cart_numerical, 0.0, dt, 0.1);
    j2_dyn.PropagateWithStm(cart_state, 0.0, dt, 0.1, stm_state);
    j2_dyn.PropagateWithStm(cart_vector, 0.0, dt, 0.1, stm_vector);

    EXPECT_NEAR_ADVEC(cart_numerical, cart_state.GetVector(), 1e-6);
    EXPECT_NEAR_ADVEC(cart_numerical, cart_vector, 1e-6);

    EXPECT_NEAR_EIGENMAT(stm_numerical, stm_state, 1e-5);
    EXPECT_NEAR_EIGENMAT(stm_numerical, stm_vector, 1e-5);
  }
}
**/
