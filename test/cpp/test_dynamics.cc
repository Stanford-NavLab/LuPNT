#include <lupnt/core/constants.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/physics/orbit_state.h>
#include <lupnt/physics/time_converter.h>

#include <catch2/catch_test_macros.hpp>
#include <iostream>

#include "utils.cc"

using namespace lupnt;

Vec6 GetClassicalOE() {
  Real a = 5740;        // [km] Semi-major axis
  Real e = 0.58;        // [-] Eccentricity
  Real i = 54.9 * DEG;  // [rad] Inclination
  Real O = 0;           // [rad] Right Ascension of Ascending Node
  Real w = 86.3;        // [rad] Argument of Perigee
  Real M = 0;           // [rad] True Anomaly

  Vec6 coe = {a, e, i, O, w, M};
  return coe;
}

MatX6 GetClassicalOEMat(int n) {
  Vec6 coe = GetClassicalOE();
  MatX6 coe_mat = coe.replicate(1, n);
  coe_mat.col(5) = VecX::LinSpaced(n, 0, TWO_PI);
  return coe_mat;
}

const double ABS_TOL = 1e-8;

TEST_CASE("TestTwoBodyDynamics") {
  // Constants
  int N_sat = 4;         // Number of satellites
  Real J2 = 0;           // [-] J2 coefficient
  Real GM = GM_MOON;     // [km^3/s^2] Gravitational parameter
  Real R_body = R_MOON;  // [km] Radius of the central body

  // Time
  int N_steps = 4;
  Real dt = 10;                                              // [s]
  Real t0 = Gregorian2Time(2024, 6, 1, 12, 45, 30);          // [s] TAI
  VecX tspan = VecX::LinSpaced(N_steps, 0, 24 * SECS_HOUR);  // [s]
  VecX tfs = t0 + tspan.array();                             // [s] TAI

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
  Ptr<IDynamics> dyn_kep_ptr = Ptr<KeplerianDynamics>(&dyn_kep);
  Ptr<IDynamics> dyn_cart_ptr = Ptr<CartesianTwoBodyDynamics>(&dyn_cart);
  Ptr<IDynamics> dyn_cart_j2_ptr = Ptr<J2CartTwoBodyDynamics>(&dyn_cart_j2);
  Ptr<IDynamics> dyn_kep_j2_ptr = Ptr<J2KeplerianDynamics>(&dyn_kep_j2);

  // Check time step
  REQUIRE_NEAR_REAL(std::static_pointer_cast<CartesianTwoBodyDynamics>(dyn_cart_ptr)->GetTimeStep(),
                    dt, ABS_TOL);
  REQUIRE_NEAR_REAL(std::static_pointer_cast<J2CartTwoBodyDynamics>(dyn_cart_j2_ptr)->GetTimeStep(),
                    dt, ABS_TOL);
  REQUIRE_NEAR_REAL(std::static_pointer_cast<J2KeplerianDynamics>(dyn_kep_j2_ptr)->GetTimeStep(),
                    dt, ABS_TOL);

  // Vec6
  Vec6 coe_vec = GetClassicalOE();
  Vec6 rv_vec = Classical2Cart(coe_vec, GM);
  Vec6 coe_j2_vec = coe_vec;
  Vec6 rv_j2_vec = rv_vec;
  Vec6 qnsoe_vec = Classical2QuasiNonsing(coe_vec, GM);
  Vec6 eqoe_vec = Classical2Equinoctial(coe_vec, GM);

  // ************************************************************************************************
  // Multple times
  // ************************************************************************************************

  // (N_steps x 6)
  MatX6 coe_prop = dyn_kep.Propagate(coe_vec, t0, tfs);
  MatX6 coe_j2_prop = dyn_kep_j2.Propagate(coe_j2_vec, t0, tfs);
  MatX6 rv_prop = dyn_cart.Propagate(rv_vec, t0, tfs);
  MatX6 rv_j2_prop = dyn_cart_j2.Propagate(rv_j2_vec, t0, tfs);
  MatX6 qnsoe_prop = dyn_kep.Propagate(qnsoe_vec, t0, tfs);
  MatX6 eqoe_prop = dyn_kep.Propagate(eqoe_vec, t0, tfs);

  // Check rows
  REQUIRE(coe_prop.rows() == N_steps);
  REQUIRE(coe_j2_prop.rows() == N_steps);
  REQUIRE(rv_prop.rows() == N_steps);
  REQUIRE(rv_j2_prop.rows() == N_steps);
  REQUIRE(qnsoe_prop.rows() == N_steps);
  REQUIRE(eqoe_prop.rows() == N_steps);

  // Check values
  MatX6 tmp;
  REQUIRE_NEAR_REAL_MAT(coe_prop, coe_j2_prop, ABS_TOL);
  tmp = Cart2Classical(rv_prop, GM);
  REQUIRE_NEAR_REAL_MAT(coe_prop, tmp, ABS_TOL);
  REQUIRE_NEAR_REAL_MAT(coe_prop, Cart2Classical(rv_j2_prop, GM), ABS_TOL);
  REQUIRE_NEAR_REAL_MAT(coe_prop, QuasiNonsing2Classical(qnsoe_prop, GM), ABS_TOL);
  REQUIRE_NEAR_REAL_MAT(coe_prop, Equinoctial2Classical(eqoe_prop, GM), ABS_TOL);

  // ************************************************************************************************
  // Multple vectors
  // ************************************************************************************************

  // (N_sat x 6)
  MatX6 coe_prop_sat = GetClassicalOEMat(N_sat);
  MatX6 coe_j2_prop_sat = coe_prop_sat;
  MatX6 rv_matx6 = Classical2Cart(coe_prop_sat, GM);
  MatX6 rv_j2_prop_sat = Classical2Cart(coe_prop_sat, GM);
  MatX6 qnsoe_matx6 = Classical2QuasiNonsing(coe_prop_sat, GM);
  MatX6 eqoe_matx6 = Classical2Equinoctial(coe_prop_sat, GM);

  // Check times
  REQUIRE_NEAR_REAL(t0, tfs(0), ABS_TOL);

  // Propagate loop
  for (int i = 1; i < N_steps; i++) {
    coe_prop_sat = dyn_kep.Propagate(coe_prop_sat, tfs(i - 1), tfs(i));
    coe_j2_prop_sat = dyn_kep_j2.Propagate(coe_j2_prop_sat, tfs(i - 1), tfs(i));
    rv_matx6 = dyn_cart.Propagate(rv_matx6, tfs(i - 1), tfs(i));
    rv_j2_prop_sat = dyn_cart_j2.Propagate(rv_j2_prop_sat, tfs(i - 1), tfs(i));
    qnsoe_matx6 = dyn_kep.Propagate(qnsoe_matx6, tfs(i - 1), tfs(i));
    eqoe_matx6 = dyn_kep.Propagate(eqoe_matx6, tfs(i - 1), tfs(i));

    // Check values
    REQUIRE_NEAR_REAL_VEC(coe_prop_sat.row(i), coe_prop.row(i), ABS_TOL);
    REQUIRE_NEAR_REAL_MAT(coe_prop_sat, coe_j2_prop_sat, ABS_TOL);
    REQUIRE_NEAR_REAL_MAT(coe_prop_sat, Cart2Classical(rv_matx6, GM), ABS_TOL);
    REQUIRE_NEAR_REAL_MAT(coe_prop_sat, Cart2Classical(rv_j2_prop_sat, GM), ABS_TOL);
    REQUIRE_NEAR_REAL_MAT(coe_prop_sat, QuasiNonsing2Classical(qnsoe_matx6, GM), ABS_TOL);
    REQUIRE_NEAR_REAL_MAT(coe_prop_sat, Equinoctial2Classical(eqoe_matx6, GM), ABS_TOL);
  }

  // ************************************************************************************************
  // Vectors and states
  // ************************************************************************************************

  // VecX
  VecX coe_vecx = coe_vec;
  VecX coe_j2_vecx = coe_j2_vec;
  VecX rv_vecx = rv_vec;
  VecX rv_j2_vecx = rv_j2_vec;
  VecX qnsoe_vecx = qnsoe_vec;
  VecX eqoe_vecx = eqoe_vec;

  // Mat6
  Mat6d coe_stm;
  Mat6d coe_j2_stm;
  Mat6d rv_stm;
  Mat6d rv_j2_stm;
  Mat6d qnsoe_stm;
  Mat6d eqoe_stm;

  // MatX
  MatXd coe_stmx(6, 6);
  MatXd coe_j2_stmx(6, 6);
  MatXd rv_stmx(6, 6);
  MatXd rv_j2_stmx(6, 6);
  MatXd qnsoe_stmx(6, 6);
  MatXd eqoe_stmx(6, 6);

  // State
  ClassicalOE coe_state(coe_vec, Frame::MOON_CI);
  ClassicalOE coe_j2_state(coe_j2_vec, Frame::MOON_CI);
  CartesianOrbitState rv_state(rv_vec, Frame::MOON_CI);
  CartesianOrbitState rv_j2_state(rv_j2_vec, Frame::MOON_CI);
  QuasiNonsingOE qnsoe_state(qnsoe_vec, Frame::MOON_CI);
  EquinoctialOE eqoe_state(eqoe_vec, Frame::MOON_CI);

  // Ptr<IState>
  Ptr<IState> coe_ptr = MakePtr<ClassicalOE>(coe_state);
  Ptr<IState> coe_j2_ptr = MakePtr<ClassicalOE>(coe_j2_state);
  Ptr<IState> rv_ptr = MakePtr<CartesianOrbitState>(rv_state);
  Ptr<IState> rv_j2_ptr = MakePtr<CartesianOrbitState>(rv_j2_state);
  Ptr<IState> qnsoe_ptr = MakePtr<QuasiNonsingOE>(qnsoe_state);
  Ptr<IState> eqoe_ptr = MakePtr<EquinoctialOE>(eqoe_state);

  // Propagate loop
  for (int i = 1; i < N_steps; i++) {
    // Vec6
    coe_vec = dyn_kep.Propagate(coe_vec, tfs(i - 1), tfs(i), &coe_stm);
    coe_j2_vec = dyn_kep_j2.Propagate(coe_j2_vec, tfs(i - 1), tfs(i), &coe_j2_stm);
    rv_vec = dyn_cart.Propagate(rv_vec, tfs(i - 1), tfs(i), &rv_stm);
    rv_j2_vec = dyn_cart_j2.Propagate(rv_j2_vec, tfs(i - 1), tfs(i), &rv_j2_stm);
    qnsoe_vec = dyn_kep.Propagate(qnsoe_vec, tfs(i - 1), tfs(i), &qnsoe_stm);
    eqoe_vec = dyn_kep.Propagate(eqoe_vec, tfs(i - 1), tfs(i), &eqoe_stm);

    REQUIRE_NEAR_REAL_VEC(coe_vec, coe_j2_vec, ABS_TOL);
    REQUIRE_NEAR_REAL_VEC(coe_vec, Cart2Classical(rv_vec, GM), ABS_TOL);
    REQUIRE_NEAR_REAL_VEC(coe_vec, Cart2Classical(rv_j2_vec, GM), ABS_TOL);
    REQUIRE_NEAR_REAL_VEC(coe_vec, QuasiNonsing2Classical(qnsoe_vec, GM), ABS_TOL);
    REQUIRE_NEAR_REAL_VEC(coe_vec, Equinoctial2Classical(eqoe_vec, GM), ABS_TOL);

    // VecX
    coe_vecx = dyn_kep.Propagate(coe_vecx, tfs(i - 1), tfs(i), &coe_stmx);
    coe_j2_vecx = dyn_kep_j2.Propagate(coe_j2_vecx, tfs(i - 1), tfs(i), &coe_j2_stmx);
    rv_vecx = dyn_cart.Propagate(rv_vecx, tfs(i - 1), tfs(i), &rv_stmx);
    rv_j2_vecx = dyn_cart_j2.Propagate(rv_j2_vecx, tfs(i - 1), tfs(i), &rv_j2_stmx);
    qnsoe_vecx = dyn_kep.Propagate(qnsoe_vecx, tfs(i - 1), tfs(i), &qnsoe_stmx);
    eqoe_vecx = dyn_kep.Propagate(eqoe_vecx, tfs(i - 1), tfs(i), &eqoe_stmx);

    REQUIRE_NEAR_REAL_VEC(coe_vecx, coe_vec, ABS_TOL);
    REQUIRE_NEAR_REAL_VEC(coe_j2_vecx, coe_j2_vec, ABS_TOL);
    REQUIRE_NEAR_REAL_VEC(rv_vecx, rv_vec, ABS_TOL);
    REQUIRE_NEAR_REAL_VEC(rv_j2_vecx, rv_j2_vec, ABS_TOL);
    REQUIRE_NEAR_REAL_VEC(qnsoe_vecx, qnsoe_vec, ABS_TOL);
    REQUIRE_NEAR_REAL_VEC(eqoe_vecx, eqoe_vec, ABS_TOL);

    REQUIRE_NEAR_DOUBLE_MAT(coe_stmx, coe_stm, ABS_TOL);
    REQUIRE_NEAR_DOUBLE_MAT(coe_j2_stmx, coe_j2_stm, ABS_TOL);
    REQUIRE_NEAR_DOUBLE_MAT(rv_stmx, rv_stm, ABS_TOL);
    REQUIRE_NEAR_DOUBLE_MAT(rv_j2_stmx, rv_j2_stm, ABS_TOL);
    REQUIRE_NEAR_DOUBLE_MAT(qnsoe_stmx, qnsoe_stm, ABS_TOL);
    REQUIRE_NEAR_DOUBLE_MAT(eqoe_stmx, eqoe_stm, ABS_TOL);

    // OrbitState
    coe_ptr = dyn_kep_ptr->PropagateState(coe_ptr, tfs(i - 1), tfs(i), &coe_stmx);
    coe_j2_ptr = dyn_kep_ptr->PropagateState(coe_j2_ptr, tfs(i - 1), tfs(i), &coe_j2_stmx);
    rv_ptr = dyn_cart_ptr->PropagateState(rv_ptr, tfs(i - 1), tfs(i), &rv_stmx);
    rv_j2_ptr = dyn_cart_j2_ptr->PropagateState(rv_j2_ptr, tfs(i - 1), tfs(i), &rv_j2_stmx);
    qnsoe_ptr = dyn_kep_ptr->PropagateState(qnsoe_ptr, tfs(i - 1), tfs(i), &qnsoe_stmx);
    eqoe_ptr = dyn_kep_ptr->PropagateState(eqoe_ptr, tfs(i - 1), tfs(i), &eqoe_stmx);

    REQUIRE_NEAR_REAL_VEC(coe_ptr->GetVec(), coe_vec, ABS_TOL);
    REQUIRE_NEAR_REAL_VEC(coe_j2_ptr->GetVec(), coe_j2_vec, ABS_TOL);
    REQUIRE_NEAR_REAL_VEC(rv_ptr->GetVec(), rv_vec, ABS_TOL);
    REQUIRE_NEAR_REAL_VEC(rv_j2_ptr->GetVec(), rv_j2_vec, ABS_TOL);
    REQUIRE_NEAR_REAL_VEC(qnsoe_ptr->GetVec(), qnsoe_vec, ABS_TOL);
    REQUIRE_NEAR_REAL_VEC(eqoe_ptr->GetVec(), eqoe_vec, ABS_TOL);

    REQUIRE_NEAR_DOUBLE_MAT(coe_stmx, coe_stm, ABS_TOL);
    REQUIRE_NEAR_DOUBLE_MAT(coe_j2_stmx, coe_j2_stm, ABS_TOL);
    REQUIRE_NEAR_DOUBLE_MAT(rv_stmx, rv_stm, ABS_TOL);
    REQUIRE_NEAR_DOUBLE_MAT(rv_j2_stmx, rv_j2_stm, ABS_TOL);
    REQUIRE_NEAR_DOUBLE_MAT(qnsoe_stmx, qnsoe_stm, ABS_TOL);
    REQUIRE_NEAR_DOUBLE_MAT(eqoe_stmx, eqoe_stm, ABS_TOL);
  }
}
