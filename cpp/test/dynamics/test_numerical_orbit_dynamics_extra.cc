#include <lupnt/core/config.h>
#include <lupnt/dynamics/numerical_orbit_dynamics.h>
#include <lupnt/environment/body.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;
using Catch::Approx;

namespace {
  // Two-body specific orbital energy and semi-major axis (SI units).
  double SpecificEnergy(const Vec6& rv, double GM) {
    double r = rv.head(3).norm().val();
    double v2 = rv.tail(3).squaredNorm().val();
    return 0.5 * v2 - GM / r;
  }
  double SemiMajorAxis(const Vec6& rv, double GM) { return -GM / (2.0 * SpecificEnergy(rv, GM)); }
}  // namespace

// ---------------------------------------------------------------------------
// ParseForceModelSpec: config-parsing branches (no SPICE runtime needed).
// ---------------------------------------------------------------------------
TEST_CASE("dynamics.numerical_orbit_dynamics_extra.parse_force_model_spec_full") {
  Config fm = YAML::Load(
      "bodies:\n"
      "  - MOON: {n: 8, m: 6}\n"
      "  - EARTH: {}\n"
      "  - SUN: {}\n"
      "relativity: true\n"
      "CR: 1.5\n"
      "area: 2.0\n"
      "mass: 100.0\n");
  ForceModelSpec spec = ParseForceModelSpec(fm);

  REQUIRE(spec.moon_degree == 8);
  REQUIRE(spec.moon_order == 6);
  REQUIRE(spec.include_earth);
  REQUIRE(spec.include_sun);
  REQUIRE(spec.relativity);
  REQUIRE(spec.has_srp);
  REQUIRE(spec.srp_cr == Approx(1.5));
  REQUIRE(spec.srp_area_m2 == Approx(2.0));
  REQUIRE(spec.srp_mass_kg == Approx(100.0));
}

TEST_CASE("dynamics.numerical_orbit_dynamics_extra.parse_force_model_spec_defaults") {
  // Empty / minimal block -> everything at defaults.
  Config fm = YAML::Load("{}");
  ForceModelSpec spec = ParseForceModelSpec(fm);
  REQUIRE(spec.moon_degree == 0);
  REQUIRE(spec.moon_order == 0);
  REQUIRE_FALSE(spec.include_earth);
  REQUIRE_FALSE(spec.include_sun);
  REQUIRE_FALSE(spec.relativity);
  REQUIRE_FALSE(spec.has_srp);
}

TEST_CASE("dynamics.numerical_orbit_dynamics_extra.parse_force_model_spec_use_relativity_alias") {
  // The `use_relativity` spelling is accepted as an alias for `relativity`.
  Config fm = YAML::Load(
      "bodies:\n"
      "  - MOON: {}\n"
      "use_relativity: true\n");
  ForceModelSpec spec = ParseForceModelSpec(fm);
  REQUIRE(spec.relativity);
  REQUIRE(spec.moon_degree == 0);  // MOON present but no n/m -> default 0
  REQUIRE_FALSE(spec.include_earth);
}

TEST_CASE("dynamics.numerical_orbit_dynamics_extra.parse_force_model_spec_invalid_body") {
  Config fm = YAML::Load(
      "bodies:\n"
      "  - NOT_A_BODY: {}\n");
  REQUIRE_THROWS(ParseForceModelSpec(fm));
}

// ---------------------------------------------------------------------------
// NBodyDynamics body management + force-switch getters/setters.
// ---------------------------------------------------------------------------
TEST_CASE("dynamics.numerical_orbit_dynamics_extra.body_management_and_switches") {
  NBodyDynamics dyn;
  dyn.SetFrame(Frame::MOON_CI);
  REQUIRE(dyn.GetFrame() == Frame::MOON_CI);
  REQUIRE(dyn.GetBodies().empty());

  dyn.AddBody(Body::Moon());
  dyn.AddBody(Body::Earth());
  REQUIRE(dyn.GetBodies().size() == 2);

  dyn.RemoveBody(Body::Earth());
  REQUIRE(dyn.GetBodies().size() == 1);
  REQUIRE(dyn.GetBodies()[0].id == BodyId::MOON);

  // Force-model switches round-trip through their getters.
  dyn.SetUseSrp(true);
  REQUIRE(dyn.GetUseSrp());
  dyn.SetUseDrag(true);
  REQUIRE(dyn.GetUseDrag());
  dyn.SetUseRelativity(false);
  REQUIRE_FALSE(dyn.GetUseRelativity());
}

// ---------------------------------------------------------------------------
// ComputeAccelerations: the decomposition sums back to ComputeRates.
// ---------------------------------------------------------------------------
TEST_CASE("dynamics.numerical_orbit_dynamics_extra.acceleration_decomposition_sums_to_rates") {
  NBodyDynamics dyn;
  dyn.SetFrame(Frame::MOON_CI);
  dyn.SetUseRelativity(true);
  dyn.AddBody(Body::Moon());
  dyn.AddBody(Body::Earth());

  Real t = 0.0;
  Cart6 state(Vec3(5.0e6, 0.0, 0.0), Vec3(0.0, 990.0, 0.0), Frame::MOON_CI);

  Vec6 rates = dyn.ComputeRates(t, state);
  Vec3 accel_from_rates = rates.tail(3);

  // Aggregate (non-decomposed) force terms.
  auto terms = dyn.ComputeAccelerations(t, state, false);
  REQUIRE_FALSE(terms.empty());
  Vec3 sum = Vec3::Zero();
  for (const auto& kv : terms) sum += kv.second;
  for (int i = 0; i < 3; i++) {
    REQUIRE_THAT(sum(i).val(), WithinAbs(accel_from_rates(i).val(), 1e-9));
  }

  // Rates position-block equals the velocity part of the state.
  for (int i = 0; i < 3; i++) {
    REQUIRE_THAT(rates(i).val(), WithinAbs(state.tail(3)(i).val(), 1e-12));
  }

  // Per-harmonic decomposition also sums to the same acceleration.
  auto terms_decomp = dyn.ComputeAccelerations(t, state, true);
  Vec3 sum_d = Vec3::Zero();
  for (const auto& kv : terms_decomp) sum_d += kv.second;
  for (int i = 0; i < 3; i++) {
    REQUIRE_THAT(sum_d(i).val(), WithinAbs(accel_from_rates(i).val(), 1e-9));
  }
}

// ---------------------------------------------------------------------------
// Propagation of a lunar orbit (Moon + Earth + Sun as third bodies).
// ---------------------------------------------------------------------------
TEST_CASE("dynamics.numerical_orbit_dynamics_extra.propagate_lunar_orbit_energy_sanity") {
  NBodyDynamics dyn;
  dyn.SetFrame(Frame::MOON_CI);
  dyn.SetUseRelativity(false);
  dyn.AddBody(Body::Moon());
  dyn.AddBody(Body::Earth());
  dyn.AddBody(Body::Sun());
  dyn.SetTimeStep(10.0);

  double GM = GM_MOON;
  double r0 = 5.0e6;               // 5000 km Moon-centered radius
  double v0 = std::sqrt(GM / r0);  // circular speed
  Cart6 x0(Vec3(r0, 0.0, 0.0), Vec3(0.0, v0, 0.0), Frame::MOON_CI);

  Real t0 = 0.0;
  Real tf = 600.0;  // short arc
  Cart6 xf = dyn.Propagate(x0, t0, tf);
  REQUIRE(xf.size() == 6);

  Vec6 rvf = xf;
  // Radius stays physically bounded and SMA is nearly conserved over the arc
  // (third-body perturbations are tiny relative to lunar central gravity).
  double rf = rvf.head(3).norm().val();
  REQUIRE(rf > 0.9 * r0);
  REQUIRE(rf < 1.1 * r0);

  double a0 = SemiMajorAxis(Vec6(x0), GM);
  double af = SemiMajorAxis(rvf, GM);
  REQUIRE(a0 == Approx(r0).epsilon(1e-6));
  REQUIRE(af == Approx(a0).epsilon(1e-3));  // < 0.1% drift
}

TEST_CASE("dynamics.numerical_orbit_dynamics_extra.propagate_ex_stm_shape_and_symplectic") {
  NBodyDynamics dyn;
  dyn.SetFrame(Frame::MOON_CI);
  dyn.SetUseRelativity(false);
  dyn.AddBody(Body::Moon());
  dyn.SetTimeStep(10.0);

  double GM = GM_MOON;
  double r0 = 5.0e6;
  double v0 = std::sqrt(GM / r0);
  Cart6 x0(Vec3(r0, 0.0, 0.0), Vec3(0.0, v0, 0.0), Frame::MOON_CI);

  MatXd stm;
  TerminationInfo info;
  State xf = dyn.PropagateExStm(x0, 0.0, 300.0, &stm, &info);

  REQUIRE(xf.size() == 6);
  REQUIRE(stm.rows() == 6);
  REQUIRE(stm.cols() == 6);
  REQUIRE_FALSE(info.terminated);
  REQUIRE(info.reason == TerminationReason::ReachedTf);

  // A conservative (Hamiltonian) flow is symplectic, so its state-transition
  // matrix has unit determinant (to integrator accuracy).
  REQUIRE(stm.determinant() == Approx(1.0).margin(1e-3));

  // PropagateEx (no STM) returns the same final state and reaches tf.
  TerminationInfo info2;
  State xf2 = dyn.PropagateEx(x0, 0.0, 300.0, &info2);
  REQUIRE(Vec6(xf2).isApprox(Vec6(xf), 1e-9));
  REQUIRE_FALSE(info2.terminated);
}

TEST_CASE("dynamics.numerical_orbit_dynamics_extra.propagate_ex_multi_epoch") {
  NBodyDynamics dyn;
  dyn.SetFrame(Frame::MOON_CI);
  dyn.SetUseRelativity(false);
  dyn.AddBody(Body::Moon());
  dyn.SetTimeStep(10.0);

  double GM = GM_MOON;
  double r0 = 5.0e6;
  double v0 = std::sqrt(GM / r0);
  Cart6 x0(Vec3(r0, 0.0, 0.0), Vec3(0.0, v0, 0.0), Frame::MOON_CI);

  VecX tfs(3);
  tfs << 0.0, 150.0, 300.0;
  TerminationInfo info;
  MatX traj = dyn.PropagateEx(x0, 0.0, tfs, &info);

  REQUIRE(traj.rows() == 3);
  REQUIRE(traj.cols() == 6);
  // First row is the initial state (tf == t0).
  REQUIRE(traj.row(0).transpose().isApprox(Vec6(x0), 1e-9));
}
