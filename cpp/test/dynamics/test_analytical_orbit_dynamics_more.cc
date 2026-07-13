#include <lupnt/dynamics/analytical_orbit_dynamics.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "../utils.cc"
#include "lupnt/numerics/math_utils.h"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // ---------------------------------------------------------------------------
  // KeplerianDynamics<ClassicalOE>: full STM structure. The closed-form STM is
  // the identity everywhere except the d(M)/d(a) sensitivity entry (5,0). The
  // existing tests only check that single entry; here we verify every other
  // entry equals the identity, and that the sensitivity has the correct sign
  // and magnitude across several arc lengths.
  // ---------------------------------------------------------------------------
  TEST_CASE("dynamics.analytical_orbit_dynamics_more.kepler_stm_structure") {
    const Real GM = GM_EARTH;
    ClassicalOE coe(Vec6(7200.0e3, 0.02, 0.3, 0.4, 0.5, 0.6), Frame::GCRF);
    KeplerianDynamics<ClassicalOE> dyn(GM);

    Real a = coe.a();
    Real n = sqrt(GM / pow(a, 3));

    for (double dt : {30.0, 600.0, 3000.0}) {
      INFO("dt = " << dt);
      MatXd stm;
      dyn.Propagate(coe, Real(0.0), Real(dt), nullptr, &stm);
      REQUIRE(stm.rows() == 6);
      REQUIRE(stm.cols() == 6);

      double expected_dMda = (-1.5 * (n / a) * dt).val();
      for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
          double expected = (i == j) ? 1.0 : 0.0;
          if (i == 5 && j == 0) expected = expected_dMda;
          INFO("entry (" << i << "," << j << ")");
          REQUIRE_THAT(stm(i, j), WithinAbs(expected, 1e-12));
        }
      }
      // d(M)/d(a) is negative for prograde motion (larger a -> slower).
      REQUIRE(expected_dMda < 0.0);
    }
  }

  // ---------------------------------------------------------------------------
  // Backward Keplerian propagation is the exact inverse of forward propagation
  // (M advances by n*dt for any signed dt).
  // ---------------------------------------------------------------------------
  TEST_CASE("dynamics.analytical_orbit_dynamics_more.kepler_forward_backward_inverse") {
    const Real GM = GM_MOON;
    ClassicalOE coe(Vec6(3000.0e3, 0.05, 0.6, 0.1, 0.2, 0.3), Frame::MOON_CI);
    KeplerianDynamics<ClassicalOE> dyn(GM);

    Real dt = 1234.0;
    State fwd = dyn.Propagate(coe, Real(0.0), dt, nullptr);
    State back = dyn.Propagate(fwd, Real(0.0), Real(-dt.val()), nullptr);

    for (int i = 0; i < 5; ++i) REQUIRE_THAT(back(i).val(), WithinAbs(coe(i).val(), 1e-6));
    REQUIRE_THAT(WrapToPi(Real(back(5) - coe(5))).val(), WithinAbs(0.0, 1e-9));
  }

  // ---------------------------------------------------------------------------
  // The parallelized multi-epoch Propagate(x0, tfs) must agree row-by-row with
  // repeated single-epoch Propagate calls, and the non-STM Propagate overload
  // must agree with the STM overload's returned state.
  // ---------------------------------------------------------------------------
  TEST_CASE("dynamics.analytical_orbit_dynamics_more.kepler_multiepoch_and_overloads") {
    const Real GM = GM_EARTH;
    ClassicalOE coe(Vec6(8000.0e3, 0.1, 0.2, 0.3, 0.4, 0.5), Frame::GCRF);
    KeplerianDynamics<ClassicalOE> dyn(GM);

    VecX tfs(4);
    tfs << 0.0, 100.0, 500.0, 2000.0;
    MatX history = dyn.Propagate(coe, tfs);
    REQUIRE(history.rows() == 4);
    REQUIRE(history.cols() == 6);

    for (int k = 0; k < tfs.size(); ++k) {
      State single = dyn.Propagate(coe, Real(0.0), Real(tfs(k)), nullptr);
      for (int i = 0; i < 6; ++i)
        REQUIRE_THAT(history(k, i).val(), WithinAbs(single(i).val(), 1e-6));
    }

    // Non-STM vs STM overload agreement.
    MatXd stm;
    State with_stm = dyn.Propagate(coe, Real(0.0), Real(500.0), nullptr, &stm);
    State without_stm = dyn.Propagate(coe, Real(0.0), Real(500.0), nullptr);
    for (int i = 0; i < 6; ++i)
      REQUIRE_THAT(without_stm(i).val(), WithinAbs(with_stm(i).val(), 1e-12));
  }

  // ---------------------------------------------------------------------------
  // ClohessyWiltshire Propagate: a zero-length interval is an identity, and the
  // returned STM for that case is the 6x6 identity. (The non-zero-interval
  // Propagate path re-solves the CW integration constants from the *reference*
  // epoch and is exercised at the ComputeMat level in the _extra suite.)
  // ---------------------------------------------------------------------------
  TEST_CASE("dynamics.analytical_orbit_dynamics_more.clohessy_wiltshire_zero_interval") {
    Real a = 7000.0e3;
    double nd = 1.1e-3;
    Real n = nd;
    ClohessyWiltshireDynamics cw(a, n);

    RelCart6 x0(Vec3(10.0, -5.0, 3.0), Vec3(0.1, -0.2, 0.05), Frame::GCRF);
    MatXd stm;
    State xf = cw.Propagate(x0, Real(5.0), Real(5.0), nullptr, &stm);
    for (int i = 0; i < 6; ++i) REQUIRE_THAT(xf(i).val(), WithinAbs(x0(i).val(), 1e-12));

    REQUIRE(stm.rows() == 6);
    REQUIRE(stm.cols() == 6);
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j) REQUIRE_THAT(stm(i, j), WithinAbs(i == j ? 1.0 : 0.0, 1e-12));
  }

  // ---------------------------------------------------------------------------
  // YamanakaAnkersen: ComputeMat and ComputeInverseMat are mutual inverses; the
  // Propagate entry point is explicitly unimplemented for dt != 0.
  // ---------------------------------------------------------------------------
  TEST_CASE("dynamics.analytical_orbit_dynamics_more.yamanaka_ankersen_mat") {
    const Real GM = GM_EARTH;
    ClassicalOE coe_c(Vec6(8000.0e3, 0.15, 0.3, 0.4, 0.5, 0.2), Frame::GCRF);
    Cart6 rv_rtn(Vec3(100.0, 50.0, -20.0), Vec3(0.1, -0.05, 0.02), Frame::GCRF);
    YamanakaAnkersenDynamics ya(coe_c, rv_rtn, GM);

    for (double t : {150.0, 900.0}) {
      INFO("t = " << t);
      MatX M = ya.ComputeMat(Real(t));
      MatX Minv = ya.ComputeInverseMat(Real(t));
      REQUIRE(M.rows() == 6);
      REQUIRE(M.cols() == 6);
      MatX prod = M * Minv;
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
          REQUIRE_THAT(prod(i, j).val(), WithinAbs(i == j ? 1.0 : 0.0, 1e-7));
    }

    RelCart6 x0(Vec3(1.0, 2.0, 3.0), Vec3(0.0, 0.0, 0.0), Frame::GCRF);
    REQUIRE_THROWS_AS(ya.Propagate(x0, Real(0.0), Real(100.0), nullptr), std::runtime_error);
  }

  // ---------------------------------------------------------------------------
  // RoeGeometricMapping: ComputeMat is functional (6x6, finite); Propagate is
  // explicitly unimplemented for dt != 0.
  // ---------------------------------------------------------------------------
  TEST_CASE("dynamics.analytical_orbit_dynamics_more.roe_geometric_mapping") {
    const Real GM = GM_EARTH;
    ClassicalOE coe_c(Vec6(7500.0e3, 0.08, 0.5, 0.6, 0.7, 0.1), Frame::GCRF);
    QuasiNonsingROE roe(Vec6(0.0, 100.0, 50.0, -30.0, 40.0, -10.0));
    RoeGeometricMappingDynamics roe_dyn(coe_c, roe, GM);

    for (double t : {0.0, 400.0, 1200.0}) {
      INFO("t = " << t);
      MatX M = roe_dyn.ComputeMat(Real(t));
      REQUIRE(M.rows() == 6);
      REQUIRE(M.cols() == 6);
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) REQUIRE(std::isfinite(M(i, j).val()));
    }

    RelCart6 x0(Vec3(1.0, 2.0, 3.0), Vec3(0.0, 0.0, 0.0), Frame::GCRF);
    REQUIRE_THROWS_AS(roe_dyn.Propagate(x0, Real(0.0), Real(100.0), nullptr), std::runtime_error);
  }
}  // namespace
