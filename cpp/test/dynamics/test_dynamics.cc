#include <lupnt/conversions/state_conversions.h>
#include <lupnt/conversions/time_conversions.h>
#include <lupnt/core/constants.h>
#include <lupnt/dynamics/analytical_orbit_dynamics.h>
#include <lupnt/dynamics/numerical_orbit_dynamics.h>

#include <catch2/catch_test_macros.hpp>

#include "../data.cc"
#include "../utils.cc"

using namespace lupnt;

TEST_CASE("dynamics.two_body") {
  const Real GM = GM_MOON;
  const Real t0 = GregorianToTime(2024, 6, 1, 12, 45, 30);
  const Real dt = 10.0;

  ClassicalOE coe0(GetClassicalOE(), Frame::MOON_CI);
  Cart6 rv0(ClassicalToCart(coe0, GM), Frame::MOON_CI);

  KeplerianDynamics dyn_kep(GM);
  CartesianTwoBodyDynamics dyn_cart(GM);
  JToCartTwoBodyDynamics dyn_cart_j2(GM, Real(0.0), Real(R_MOON));

  dyn_cart.SetTimeStep(dt);
  dyn_cart_j2.SetTimeStep(dt);
  RequireNear(dyn_cart.GetTimeStep(), dt, ABS_TOL);
  RequireNear(dyn_cart_j2.GetTimeStep(), dt, ABS_TOL);

  SECTION("single-epoch propagation preserves compatible Keplerian and Cartesian states") {
    const Real tf = t0 + 60.0;

    ClassicalOE coe = dyn_kep.Propagate(coe0, t0, tf);
    Cart6 rv = dyn_cart.Propagate(rv0, t0, tf);
    Cart6 rv_j2 = dyn_cart_j2.Propagate(rv0, t0, tf);

    RequireNear(coe, CartToClassical(rv, GM), 1e-3);
    RequireNear(rv, rv_j2, 1e-3);
  }

  SECTION("multi-time propagation returns one row per requested epoch") {
    VecX tfs(4);
    tfs << t0, t0 + 60.0, t0 + 120.0, t0 + 180.0;

    MatX coe_prop = dyn_kep.Propagate(coe0, tfs);
    MatX rv_prop = dyn_cart.Propagate(rv0, tfs);

    REQUIRE(coe_prop.rows() == tfs.size());
    REQUIRE(coe_prop.cols() == coe0.size());
    REQUIRE(rv_prop.rows() == tfs.size());
    REQUIRE(rv_prop.cols() == rv0.size());

    RequireNear(coe_prop.row(0).transpose(), coe0, ABS_TOL);
    RequireNear(rv_prop.row(0).transpose(), rv0, ABS_TOL);
  }

  SECTION("STM overload fills a state transition matrix") {
    MatXd stm(6, 6);
    State coe = dyn_kep.Propagate(coe0, t0, t0 + 60.0, nullptr, &stm);

    REQUIRE(coe.size() == 6);
    REQUIRE(stm.rows() == 6);
    REQUIRE(stm.cols() == 6);
  }
}
