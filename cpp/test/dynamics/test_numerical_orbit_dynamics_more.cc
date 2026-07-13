#include <lupnt/core/config.h>
#include <lupnt/core/constants.h>
#include <lupnt/dynamics/numerical_orbit_dynamics.h>
#include <lupnt/environment/body.h>
#include <lupnt/numerics/integrator.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;
using Catch::Approx;

namespace {
  double SpecificEnergy(const Vec6& rv, double GM) {
    double r = rv.head(3).norm().val();
    double v2 = rv.tail(3).squaredNorm().val();
    return 0.5 * v2 - GM / r;
  }
}  // namespace

// ---------------------------------------------------------------------------
// CartesianTwoBodyDynamics: an isolated point-mass integrator. Covers the
// constructor, ComputeRates, and the inherited fixed-step Propagate. A circular
// orbit conserves both radius and specific orbital energy over the arc.
// ---------------------------------------------------------------------------
TEST_CASE("dynamics.numerical_orbit_dynamics_more.cartesian_two_body_circular_energy") {
  double GM = GM_EARTH;
  CartesianTwoBodyDynamics dyn(GM);
  dyn.SetTimeStep(5.0);

  double r0 = 7.0e6;               // 7000 km geocentric radius
  double v0 = std::sqrt(GM / r0);  // circular speed
  Cart6 x0(Vec3(r0, 0.0, 0.0), Vec3(0.0, v0, 0.0), Frame::GCRF);

  // ComputeRates: position-block returns velocity, velocity-block points inward.
  Vec6 rates = dyn.ComputeRates(0.0, x0);
  for (int i = 0; i < 3; i++) REQUIRE_THAT(rates(i).val(), WithinAbs(x0.tail(3)(i).val(), 1e-12));
  REQUIRE(rates(3).val() < 0.0);  // radial acceleration toward the center (-x)
  REQUIRE_THAT(rates(4).val(), WithinAbs(0.0, 1e-6));

  Cart6 xf = dyn.Propagate(x0, 0.0, 600.0);
  Vec6 rvf = xf;
  double rf = rvf.head(3).norm().val();
  REQUIRE_THAT(rf, WithinRel(r0, 1e-4));

  double e0 = SpecificEnergy(Vec6(x0), GM);
  double ef = SpecificEnergy(rvf, GM);
  REQUIRE(ef == Approx(e0).epsilon(1e-6));

  // Zero-length propagation is an identity (abs(tf - t0) < EPS short-circuit).
  Cart6 xsame = dyn.Propagate(x0, 10.0, 10.0);
  REQUIRE(Vec6(xsame).isApprox(Vec6(x0), 1e-12));
}

// ---------------------------------------------------------------------------
// The concrete integrators selectable through SetIntegrator all integrate the
// same two-body flow to (approximately) the same circular-orbit final state.
// This exercises the RK8 / RKF45 / PD45 branches of SetIntegrator.
// ---------------------------------------------------------------------------
TEST_CASE("dynamics.numerical_orbit_dynamics_more.set_integrator_variants") {
  double GM = GM_EARTH;
  double r0 = 7.0e6;
  double v0 = std::sqrt(GM / r0);
  Cart6 x0(Vec3(r0, 0.0, 0.0), Vec3(0.0, v0, 0.0), Frame::GCRF);

  for (IntegratorType integ :
       {IntegratorType::RK4, IntegratorType::RK8, IntegratorType::RKF45, IntegratorType::PD45}) {
    CartesianTwoBodyDynamics dyn(GM);
    dyn.SetTimeStep(5.0);
    dyn.SetIntegrator(integ);
    REQUIRE(dyn.GetIntegrator() != nullptr);
    Cart6 xf = dyn.Propagate(x0, 0.0, 300.0);
    Vec6 rvf = xf;
    for (int i = 0; i < 6; i++) REQUIRE(std::isfinite(rvf(i).val()));
    REQUIRE_THAT(rvf.head(3).norm().val(), WithinRel(r0, 1e-3));
  }
}

// ---------------------------------------------------------------------------
// JToCartTwoBodyDynamics: the J2 == 0 shortcut returns pure two-body rates
// without touching any frame machinery; the J2 != 0 path adds an oblateness
// perturbation (SPICE frame rotation) that leaves the position-block untouched.
// ---------------------------------------------------------------------------
TEST_CASE("dynamics.numerical_orbit_dynamics_more.j2cart_two_body_rates") {
  double GM = GM_EARTH;
  double R = R_EARTH;
  double J2 = 1.08262668e-3;

  Cart6 x(Vec3(7.0e6, 1.0e6, 2.0e6), Vec3(0.0, 7.0e3, 1.0e3), Frame::GCRF);

  SECTION("J2 == 0 yields the plain point-mass acceleration") {
    JToCartTwoBodyDynamics dyn(GM, 0.0, R);
    Vec6 rates = dyn.ComputeRates(0.0, x);
    Vec3 r = x.head(3);
    Vec3 a_expected = -GM * r / std::pow(r.norm().val(), 3);
    for (int i = 0; i < 3; i++) {
      REQUIRE_THAT(rates(i).val(), WithinAbs(x.tail(3)(i).val(), 1e-12));
      REQUIRE_THAT(rates(3 + i).val(), WithinAbs(a_expected(i).val(), 1e-9));
    }
  }

  SECTION("J2 != 0 perturbs the acceleration but not the position rates") {
    JToCartTwoBodyDynamics dyn(GM, J2, R, Frame::GCRF, Frame::ITRF);
    REQUIRE(dyn.GetFrame() == Frame::GCRF);
    REQUIRE(dyn.GetBodyFixedFrame() == Frame::ITRF);

    Vec6 rates = dyn.ComputeRates(0.0, x);
    for (int i = 0; i < 3; i++) REQUIRE_THAT(rates(i).val(), WithinAbs(x.tail(3)(i).val(), 1e-12));

    // The J2 term is a small correction on top of central gravity: the total
    // acceleration must stay finite and remain close to the point-mass value.
    Vec3 r = x.head(3);
    Vec3 a_pm = -GM * r / std::pow(r.norm().val(), 3);
    for (int i = 0; i < 3; i++) {
      REQUIRE(std::isfinite(rates(3 + i).val()));
      REQUIRE(std::abs((rates(3 + i) - a_pm(i)).val()) < 0.02 * a_pm.norm().val() + 1e-6);
    }
  }
}

// ---------------------------------------------------------------------------
// J2KeplerianDynamics: mean-element (COE) secular rates. Semi-major axis,
// eccentricity, and inclination have zero secular drift; the node regresses for
// a prograde orbit and the mean anomaly advances.
// ---------------------------------------------------------------------------
TEST_CASE("dynamics.numerical_orbit_dynamics_more.j2_keplerian_secular_rates") {
  double GM = GM_EARTH;
  double R = R_EARTH;
  double J2 = 1.08262668e-3;
  J2KeplerianDynamics dyn(GM, J2, R);

  double a = 7.0e6, e = 0.01, inc = 0.9, Omega = 0.3, w = 0.5, M = 0.1;
  ClassicalOE coe({a, e, inc, Omega, w, M}, Frame::GCRF);

  Vec6 rates = dyn.ComputeRates(0.0, coe);
  REQUIRE_THAT(rates(0).val(), WithinAbs(0.0, 1e-12));  // da/dt
  REQUIRE_THAT(rates(1).val(), WithinAbs(0.0, 1e-12));  // de/dt
  REQUIRE_THAT(rates(2).val(), WithinAbs(0.0, 1e-12));  // di/dt
  REQUIRE(rates(3).val() < 0.0);                        // RAAN regresses (cos i > 0)
  REQUIRE(rates(5).val() > 0.0);                        // mean anomaly advances
  for (int i = 0; i < 6; i++) REQUIRE(std::isfinite(rates(i).val()));

  // A retrograde orbit reverses the sign of the nodal rate (cos i < 0).
  ClassicalOE coe_retro({a, e, 2.6, Omega, w, M}, Frame::GCRF);
  Vec6 rates_retro = dyn.ComputeRates(0.0, coe_retro);
  REQUIRE(rates_retro(3).val() > 0.0);
}

// ---------------------------------------------------------------------------
// MoonMeanDynamics: lunar mean-element secular rates driven by J2 and the
// Earth third-body term. Semi-major axis is unchanged; every rate is finite.
// ---------------------------------------------------------------------------
TEST_CASE("dynamics.numerical_orbit_dynamics_more.moon_mean_rates") {
  MoonMeanDynamics dyn;
  double a = 6.0e6, e = 0.05, inc = 0.8, Omega = 0.2, w = 0.4, M = 0.0;
  ClassicalOE coe({a, e, inc, Omega, w, M}, Frame::MOON_CI);

  Vec6 rates = dyn.ComputeRates(0.0, coe);
  REQUIRE_THAT(rates(0).val(), WithinAbs(0.0, 1e-12));  // da/dt == 0 by construction
  for (int i = 0; i < 6; i++) REQUIRE(std::isfinite(rates(i).val()));
  // The mean longitude rate is dominated by the mean motion and stays positive.
  REQUIRE(rates(5).val() > 0.0);
}

// ---------------------------------------------------------------------------
// NBodyDynamics constructed from a unified force-model YAML config. This drives
// the config constructor: frame, integrator/dt, SRP (CR/area/mass), bodies with
// gravity degree/order, and the relativity/autodiff switches.
// ---------------------------------------------------------------------------
TEST_CASE("dynamics.numerical_orbit_dynamics_more.nbody_config_constructor_srp") {
  Config cfg = YAML::Load(
      "frame: MOON_CI\n"
      "integrator: RK4\n"
      "dt: 30.0\n"
      "relativity: false\n"
      "autodiff: false\n"
      "area: 2.0\n"
      "mass: 100.0\n"
      "CR: 1.5\n"
      "bodies:\n"
      "  - MOON: {n: 4, m: 4}\n"
      "  - EARTH: {}\n");
  NBodyDynamics dyn(cfg);

  REQUIRE(dyn.GetFrame() == Frame::MOON_CI);
  REQUIRE(dyn.GetBodies().size() == 2);
  REQUIRE(dyn.GetUseSrp());
  REQUIRE_FALSE(dyn.GetUseDrag());
  REQUIRE_FALSE(dyn.GetUseRelativity());
  REQUIRE(dyn.GetTimeStep().val() == Approx(30.0));
  // bcoeff_srp = CR * area / mass.
  REQUIRE(dyn.GetSrpCoeff().val() == Approx(1.5 * 2.0 / 100.0));
  REQUIRE(dyn.GetUnits().length == Approx(SI_UNITS.length));
}

TEST_CASE("dynamics.numerical_orbit_dynamics_more.nbody_config_constructor_km_drag") {
  Config cfg = YAML::Load(
      "frame: GCRF\n"
      "integrator: RKF45\n"
      "abstol: 1.0e-9\n"
      "reltol: 1.0e-9\n"
      "max_iter: 40\n"
      "units: km_s_kg\n"
      "area: 3.0\n"
      "mass: 200.0\n"
      "CD: 2.2\n"
      "bodies:\n"
      "  - EARTH: {n: 2, m: 0}\n");
  NBodyDynamics dyn(cfg);

  REQUIRE(dyn.GetFrame() == Frame::GCRF);
  REQUIRE(dyn.GetUseDrag());
  REQUIRE_FALSE(dyn.GetUseSrp());
  REQUIRE(dyn.GetDragCoeff().val() == Approx(2.2 * 3.0 / 200.0));
  REQUIRE(dyn.GetUnits().length == Approx(KILOMETER));
  REQUIRE(dyn.GetBodies().size() == 1);
}

// ---------------------------------------------------------------------------
// Direct SRP/drag coefficient setters and their getters, plus the guards on
// SetUnits and the autodiff-gated STM Propagate overload.
// ---------------------------------------------------------------------------
TEST_CASE("dynamics.numerical_orbit_dynamics_more.coefficient_setters_and_guards") {
  NBodyDynamics dyn;
  dyn.SetFrame(Frame::MOON_CI);

  dyn.SetSrpCoefficient(1.8, 4.0, 50.0);
  REQUIRE(dyn.GetUseSrp());
  REQUIRE(dyn.GetSrpCoeff().val() == Approx(1.8 * 4.0 / 50.0));

  dyn.SetDragCoefficient(2.0, 5.0, 25.0);
  REQUIRE(dyn.GetUseDrag());
  REQUIRE(dyn.GetDragCoeff().val() == Approx(2.0 * 5.0 / 25.0));

  dyn.SetSrpCoeff(0.05);
  REQUIRE(dyn.GetSrpCoeff().val() == Approx(0.05));
  dyn.SetDragCoeff(0.07);
  REQUIRE(dyn.GetDragCoeff().val() == Approx(0.07));

  // SetUnits must be called before any body is added.
  dyn.AddBody(Body::Moon());
  REQUIRE_THROWS(dyn.SetUnits(KM_S_KG_UNITS));

  // The STM-producing Propagate overload requires autodiff to be enabled.
  Cart6 x0(Vec3(5.0e6, 0.0, 0.0), Vec3(0.0, 990.0, 0.0), Frame::MOON_CI);
  MatXd stm;
  REQUIRE(dyn.GetUseRelativity());  // default relativity is on
  REQUIRE_THROWS(dyn.Propagate(x0, 0.0, 100.0, nullptr, &stm));
}

// ---------------------------------------------------------------------------
// ComputeAccelerations exposes SRP and drag as their own labelled terms when
// those forces are enabled (Earth + Moon so drag has an atmosphere to act on).
// ---------------------------------------------------------------------------
TEST_CASE("dynamics.numerical_orbit_dynamics_more.acceleration_terms_include_srp_and_drag") {
  NBodyDynamics dyn;
  dyn.SetFrame(Frame::GCRF);
  dyn.SetUseRelativity(false);
  dyn.AddBody(Body::Earth());
  dyn.SetSrpCoefficient(1.5, 2.0, 100.0);
  dyn.SetDragCoefficient(2.2, 2.0, 100.0);

  // A low Earth orbit so the drag model has a non-trivial atmosphere.
  double r0 = R_EARTH + 3.0e5;  // ~300 km altitude
  double v0 = std::sqrt(GM_EARTH / r0);
  Cart6 x(Vec3(r0, 0.0, 0.0), Vec3(0.0, v0, 0.0), Frame::GCRF);

  auto terms = dyn.ComputeAccelerations(0.0, x, false);
  REQUIRE(terms.count("srp") == 1);
  REQUIRE(terms.count("drag") == 1);
  REQUIRE(terms.count("EARTH_gravity") == 1);
  for (const auto& kv : terms) {
    for (int i = 0; i < 3; i++) REQUIRE(std::isfinite(kv.second(i).val()));
  }
  // Both perturbations contribute a non-zero acceleration.
  REQUIRE(terms.at("drag").norm().val() > 0.0);
  REQUIRE(terms.at("srp").norm().val() > 0.0);
}

// ---------------------------------------------------------------------------
// RemoveBody on an absent body is a no-op, and AddBody rejects a duplicate.
// ---------------------------------------------------------------------------
TEST_CASE("dynamics.numerical_orbit_dynamics_more.body_add_remove_edge_cases") {
  NBodyDynamics dyn;
  dyn.SetFrame(Frame::MOON_CI);
  dyn.AddBody(Body::Moon());

  // Removing a body that was never added leaves the set unchanged.
  dyn.RemoveBody(Body::Earth());
  REQUIRE(dyn.GetBodies().size() == 1);

  // Adding the same body twice is rejected.
  REQUIRE_THROWS(dyn.AddBody(Body::Moon()));

  dyn.RemoveBody(Body::Moon());
  REQUIRE(dyn.GetBodies().empty());
}
