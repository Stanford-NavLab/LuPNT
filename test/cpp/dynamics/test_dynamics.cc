#include <lupnt/core/constants.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/physics/orbit_state.h>
#include <lupnt/physics/time_converter.h>

#include <catch2/catch_test_macros.hpp>
#include <iostream>

#include "../utils.cc"

using namespace lupnt;

Vec6 GetClassicalOE() {
  Real a = 5740;        // [km] Semi-major axis
  Real e = 0.58;        // [-] Eccentricity
  Real i = 54.9 * RAD;  // [rad] Inclination
  Real O = 0;           // [rad] Right Ascension of Ascending Node
  Real w = 86.3 * RAD;  // [rad] Argument of Perigee
  Real M = 0;           // [rad] True Anomaly

  Vec6 coe = {a, e, i, O, w, M};
  return coe;
}

MatX6 GetClassicalOEMat(int n) {
  Vec6 coe = GetClassicalOE();
  MatX6 coe_mat = coe.transpose().replicate(n, 1);
  coe_mat.col(5) = VecX::LinSpaced(n, 0, TWO_PI);
  return coe_mat;
}

const double ABS_TOL = 1e-6;
const double REL_TOL = 1e-6;

TEST_CASE("TestTwoBodyDynamics") {
  // Constants
  int N_sat = 3;         // Number of satellites
  Real J2 = 0;           // [-] J2 coefficient
  Real GM = GM_MOON;     // [km^3/s^2] Gravitational parameter
  Real R_body = R_MOON;  // [km] Radius of the central body

  // Time
  int N_steps = 4;
  Real dt = 10;                                      // [s] Integration time step
  Real Dt = 6 * SECS_HOUR;                           // [s] Total integration time
  Real t0 = Gregorian2Time(2024, 6, 1, 12, 45, 30);  // [s] Start time, TAI
  VecX tspan = VecX::LinSpaced(N_steps, 0, Dt);      // [s] Time span
  VecX tfs = t0 + tspan.array();                     // [s] Times, TAI

  // Dynamics
  KeplerianDynamics dyn_kep(GM);
  CartesianTwoBodyDynamics dyn_cart(GM, IntegratorType::RK4);
  J2CartTwoBodyDynamics dyn_cart_j2(GM, J2, R_body, IntegratorType::RK4);
  J2KeplerianDynamics dyn_kep_j2(GM, J2, R_body, IntegratorType::RK4);

  // Time step
  dyn_cart.SetTimeStep(dt);
  dyn_cart_j2.SetTimeStep(dt);
  dyn_kep_j2.SetTimeStep(dt);

  // IDynamics
  Ptr<IDynamics> dyn_kep_ptr = MakePtr<KeplerianDynamics>(dyn_kep);
  Ptr<IDynamics> dyn_cart_ptr = MakePtr<CartesianTwoBodyDynamics>(dyn_cart);
  Ptr<IDynamics> dyn_cart_j2_ptr = MakePtr<J2CartTwoBodyDynamics>(dyn_cart_j2);
  Ptr<IDynamics> dyn_kep_j2_ptr = MakePtr<J2KeplerianDynamics>(dyn_kep_j2);

  // Check time step
  RequireNear(std::static_pointer_cast<CartesianTwoBodyDynamics>(dyn_cart_ptr)->GetTimeStep(), dt,
              ABS_TOL);
  RequireNear(std::static_pointer_cast<J2CartTwoBodyDynamics>(dyn_cart_j2_ptr)->GetTimeStep(), dt,
              ABS_TOL);
  RequireNear(std::static_pointer_cast<J2KeplerianDynamics>(dyn_kep_j2_ptr)->GetTimeStep(), dt,
              ABS_TOL);

  // Vec6
  Vec6 coe_vec = GetClassicalOE();
  Vec6 rv_vec = Classical2Cart(coe_vec, GM);
  Vec6 coe_j2_vec = coe_vec;
  Vec6 rv_j2_vec = rv_vec;

  // ************************************************************************************************
  // Multple times
  // ************************************************************************************************

  // (N_steps x 6)
  MatX6 coe_prop = dyn_kep.Propagate(coe_vec, t0, tfs);
  MatX6 coe_j2_prop = dyn_kep_j2.Propagate(coe_j2_vec, t0, tfs);
  MatX6 rv_prop = dyn_cart.Propagate(rv_vec, t0, tfs);
  MatX6 rv_j2_prop = dyn_cart_j2.Propagate(rv_j2_vec, t0, tfs);

  coe_j2_prop.col(5) = Wrap2Pi(coe_j2_prop.col(5));

  // Check rows
  REQUIRE(coe_prop.rows() == N_steps);
  REQUIRE(coe_j2_prop.rows() == N_steps);
  REQUIRE(rv_prop.rows() == N_steps);
  REQUIRE(rv_j2_prop.rows() == N_steps);

  // Check values
  RequireNear(coe_prop, coe_j2_prop, ABS_TOL);
  RequireNear(coe_prop, Cart2Classical(rv_prop, GM), ABS_TOL);
  RequireNear(coe_prop, Cart2Classical(rv_j2_prop, GM), ABS_TOL);

  // ************************************************************************************************
  // Multple vectors
  // ************************************************************************************************

  // (N_sat x 6)
  MatX6 coe_prop_sat = GetClassicalOEMat(N_sat);
  MatX6 coe_j2_prop_sat = coe_prop_sat;
  MatX6 rv_matx6 = Classical2Cart(coe_prop_sat, GM);
  MatX6 rv_j2_prop_sat = Classical2Cart(coe_prop_sat, GM);

  // Check times
  RequireNear(t0, tfs(0), ABS_TOL);

  // Propagate loop
  for (int i = 1; i < N_steps; i++) {
    coe_prop_sat = dyn_kep.Propagate(coe_prop_sat, tfs(i - 1), tfs(i));
    coe_j2_prop_sat = dyn_kep_j2.Propagate(coe_j2_prop_sat, tfs(i - 1), tfs(i));
    rv_matx6 = dyn_cart.Propagate(rv_matx6, tfs(i - 1), tfs(i));
    rv_j2_prop_sat = dyn_cart_j2.Propagate(rv_j2_prop_sat, tfs(i - 1), tfs(i));

    coe_j2_prop_sat.col(5) = Wrap2Pi(coe_j2_prop_sat.col(5));

    // Check values
    RequireNear(coe_prop_sat.row(0), coe_prop.row(i), ABS_TOL);
    RequireNear(coe_prop_sat, coe_j2_prop_sat, ABS_TOL);
    RequireNear(coe_prop_sat, Cart2Classical(rv_matx6, GM), ABS_TOL);
    RequireNear(coe_prop_sat, Cart2Classical(rv_j2_prop_sat, GM), ABS_TOL);
  }

  // ************************************************************************************************
  // Vectors and states
  // ************************************************************************************************

  // VecX
  VecX coe_vecx = coe_vec;
  VecX coe_j2_vecx = coe_j2_vec;
  VecX rv_vecx = rv_vec;
  VecX rv_j2_vecx = rv_j2_vec;

  // Mat6
  Mat6d coe_stm;
  Mat6d coe_j2_stm;
  Mat6d rv_stm;
  Mat6d rv_j2_stm;

  // MatX
  MatXd coe_stmx(6, 6);
  MatXd coe_j2_stmx(6, 6);
  MatXd rv_stmx(6, 6);
  MatXd rv_j2_stmx(6, 6);

  // State
  ClassicalOE coe_state(coe_vec, Frame::MOON_CI);
  ClassicalOE coe_j2_state(coe_j2_vec, Frame::MOON_CI);
  CartesianOrbitState rv_state(rv_vec, Frame::MOON_CI);
  CartesianOrbitState rv_j2_state(rv_j2_vec, Frame::MOON_CI);

  // Ptr<IState>
  Ptr<IState> coe_ptr = MakePtr<ClassicalOE>(coe_state);
  Ptr<IState> coe_j2_ptr = MakePtr<ClassicalOE>(coe_j2_state);
  Ptr<IState> rv_ptr = MakePtr<CartesianOrbitState>(rv_state);
  Ptr<IState> rv_j2_ptr = MakePtr<CartesianOrbitState>(rv_j2_state);

  // Propagate loop
  for (int i = 1; i < N_steps; i++) {
    // Vec6
    coe_vec = dyn_kep.Propagate(coe_vec, tfs(i - 1), tfs(i), &coe_stm);
    coe_j2_vec = dyn_kep_j2.Propagate(coe_j2_vec, tfs(i - 1), tfs(i), &coe_j2_stm);
    rv_vec = dyn_cart.Propagate(rv_vec, tfs(i - 1), tfs(i), &rv_stm);
    rv_j2_vec = dyn_cart_j2.Propagate(rv_j2_vec, tfs(i - 1), tfs(i), &rv_j2_stm);

    coe_j2_vec(5) = Wrap2Pi(coe_j2_vec(5));

    RequireNear(coe_vec, coe_j2_vec, ABS_TOL);
    RequireNear(coe_vec, Cart2Classical(rv_vec, GM), ABS_TOL);
    RequireNear(coe_vec, Cart2Classical(rv_j2_vec, GM), ABS_TOL);

    // VecX
    coe_vecx = dyn_kep.Propagate(coe_vecx, tfs(i - 1), tfs(i), &coe_stmx);
    coe_j2_vecx = dyn_kep_j2.Propagate(coe_j2_vecx, tfs(i - 1), tfs(i), &coe_j2_stmx);
    rv_vecx = dyn_cart.Propagate(rv_vecx, tfs(i - 1), tfs(i), &rv_stmx);
    rv_j2_vecx = dyn_cart_j2.Propagate(rv_j2_vecx, tfs(i - 1), tfs(i), &rv_j2_stmx);

    coe_j2_vecx(5) = Wrap2Pi(coe_j2_vecx(5));

    RequireNear(coe_vecx, coe_vec, ABS_TOL);
    RequireNear(coe_j2_vecx, coe_j2_vec, ABS_TOL);
    RequireNear(rv_vecx, rv_vec, ABS_TOL);
    RequireNear(rv_j2_vecx, rv_j2_vec, ABS_TOL);

    RequireNear(coe_stmx, coe_stm, ABS_TOL);
    RequireNear(coe_j2_stmx, coe_j2_stm, ABS_TOL);
    RequireNear(rv_stmx, rv_stm, ABS_TOL);
    RequireNear(rv_j2_stmx, rv_j2_stm, ABS_TOL);

    // OrbitState
    coe_ptr = dyn_kep_ptr->PropagateState(coe_ptr, tfs(i - 1), tfs(i), &coe_stmx);
    coe_j2_ptr = dyn_kep_ptr->PropagateState(coe_j2_ptr, tfs(i - 1), tfs(i), &coe_j2_stmx);
    rv_ptr = dyn_cart_ptr->PropagateState(rv_ptr, tfs(i - 1), tfs(i), &rv_stmx);
    rv_j2_ptr = dyn_cart_j2_ptr->PropagateState(rv_j2_ptr, tfs(i - 1), tfs(i), &rv_j2_stmx);

    coe_j2_ptr->SetValue(5, Wrap2Pi(coe_j2_ptr->GetValue(5)));

    RequireNear(coe_ptr->GetVec(), coe_vec, ABS_TOL);
    RequireNear(coe_j2_ptr->GetVec(), coe_j2_vec, ABS_TOL);
    RequireNear(rv_ptr->GetVec(), rv_vec, ABS_TOL);
    RequireNear(rv_j2_ptr->GetVec(), rv_j2_vec, ABS_TOL);

    RequireNear(coe_stmx, coe_stm, ABS_TOL);
    RequireNear(coe_j2_stmx, coe_j2_stm, ABS_TOL);
    RequireNear(rv_stmx, rv_stm, ABS_TOL);
    RequireNear(rv_j2_stmx, rv_j2_stm, ABS_TOL);
  }
}
