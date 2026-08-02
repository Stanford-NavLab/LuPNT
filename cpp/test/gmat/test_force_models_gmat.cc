#include <lupnt/conversions/epoch.h>
#include <lupnt/core/constants.h>
#include <lupnt/dynamics/numerical_orbit_dynamics.h>
#include <lupnt/environment/body.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// ============================================================================
// Cross-validation of force-model PROPAGATION against GMAT R2026a.
//
// GMAT reports states, not accelerations, so these are propagation-level
// comparisons (like `dynamics.j2_propagation_gmat_reference`): a spacecraft is
// propagated with one force model active and the Cartesian state after a fixed
// arc is compared. The reference is `gmat_forces_reference.json`, produced by
// `gen_gmat_forces_reference.py` against GMAT R2026a pointed at LuPNT's DE440,
// with LuPNT's constants and the shared EGM96 `.cof` forced onto GMAT.
//
// LuPNT's individual force-model *algorithms* are additionally checked at the
// formula level against Orekit (see test_force_models_orekit.cc); these GMAT
// tests exercise the assembled propagation (force model + integrator + frame).
// ============================================================================

namespace {
  Vec6 ReadVec6G(const nlohmann::json& arr) {
    Vec6 v;
    for (int i = 0; i < 6; i++) v(i) = arr.at(i).get<double>();
    return v;
  }

  // 2024-03-15T12:00:00 UTC as TDB seconds past J2000 (the fixture epoch).
  double ScenarioEpochTdb() {
    Epoch e = Epoch::FromGregorian(2024, 3, 15, 12, 0, 0.0, Time::UTC);
    return e.To(Time::TDB).ToSeconds().val();
  }

  // Propagate a Cartesian GCRF state from the scenario epoch by `dt` seconds.
  Vec6 Propagate(NBodyDynamics& dyn, const Vec6& cart0_v, double t0, double dt) {
    Cart6 cart0(cart0_v, Frame::GCRF);
    dyn.SetIntegrator(IntegratorType::RK8);
    dyn.SetTimeStep(10.0);
    State xf = dyn.Propagate(cart0, Real(t0), Real(t0 + dt));
    Cart6 cf(xf);
    Vec6 out;
    for (int i = 0; i < 6; i++) out(i) = cf(i);
    return out;
  }
}  // namespace

// ----------------------------------------------------------------------------
// 8x8 spherical-harmonic gravity (EGM96, shared .cof on both sides)
//
// Both sides use LuPNT's EGM96.cof, DE440, and LuPNT's GM/R, so a residual is
// integrator + Earth-orientation (ITRF vs GMAT EarthFixed) difference only.
// Observed ~0.19 m over 2 h.
// ----------------------------------------------------------------------------
TEST_CASE("dynamics.gravity_propagation_gmat_reference") {
  nlohmann::json data = LoadTestJson("gmat/data/gmat_forces_reference.json");
  const auto& sec = data["gravity_propagation"];
  const int n = sec["n_max"].get<int>();
  const int m = sec["m_max"].get<int>();
  const double t0 = ScenarioEpochTdb();
  constexpr double pos_tol = 3.0;   // [m]
  constexpr double vel_tol = 5e-3;  // [m/s]

  for (const auto& kase : sec["cases"]) {
    Vec6 cart0 = ReadVec6G(kase["cart0"]);
    for (const auto& pc : kase["propagated"]) {
      double dt = pc["dt"].get<double>();
      DYNAMIC_SECTION("dt = " << dt << " s") {
        NBodyDynamics dyn;
        dyn.SetFrame(Frame::GCRF);
        dyn.AddBody(Body::Earth(n, m));
        dyn.SetUseRelativity(false);
        Vec6 xf = Propagate(dyn, cart0, t0, dt);
        Vec6 ref = ReadVec6G(pc["cart"]);
        for (int i = 0; i < 3; i++) {
          RequireNear(xf(i), ref(i), pos_tol);
          RequireNear(xf(i + 3), ref(i + 3), vel_tol);
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------
// Third-body point-mass perturbation (Sun + Moon, DE440 on both sides)
//
// GMAT PointMasses = {Earth, Sun, Luna}; LuPNT adds Earth (central point mass),
// Sun and Moon as third bodies. Same DE440, so ~0.004 m over 2 h.
// ----------------------------------------------------------------------------
TEST_CASE("dynamics.third_body_propagation_gmat_reference") {
  nlohmann::json data = LoadTestJson("gmat/data/gmat_forces_reference.json");
  const auto& sec = data["third_body_propagation"];
  const double t0 = ScenarioEpochTdb();
  constexpr double pos_tol = 0.2;   // [m]
  constexpr double vel_tol = 1e-3;  // [m/s]

  for (const auto& kase : sec["cases"]) {
    Vec6 cart0 = ReadVec6G(kase["cart0"]);
    for (const auto& pc : kase["propagated"]) {
      double dt = pc["dt"].get<double>();
      DYNAMIC_SECTION("dt = " << dt << " s") {
        NBodyDynamics dyn;
        dyn.SetFrame(Frame::GCRF);
        dyn.AddBody(Body::Earth());  // central point mass
        dyn.AddBody(Body::Sun());
        dyn.AddBody(Body::Moon());
        dyn.SetUseRelativity(false);
        Vec6 xf = Propagate(dyn, cart0, t0, dt);
        Vec6 ref = ReadVec6G(pc["cart"]);
        for (int i = 0; i < 3; i++) {
          RequireNear(xf(i), ref(i), pos_tol);
          RequireNear(xf(i + 3), ref(i + 3), vel_tol);
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------
// Solar radiation pressure (cannonball), high sunward orbit (no eclipse)
//
// LuPNT's SRP takes the Sun position from SPICE internally, so the Sun is NOT
// added as a gravitating body (which would add third-body gravity GMAT's
// SrpFM = PointMasses{Earth} + SRP does not have). Matched flux (1361 W/m^2),
// Cr, area, mass. The orbit is chosen sunlit so the two shadow models never
// engage. Observed ~0.003 m over 2 h.
// ----------------------------------------------------------------------------
TEST_CASE("dynamics.srp_propagation_gmat_reference") {
  nlohmann::json data = LoadTestJson("gmat/data/gmat_forces_reference.json");
  const auto& sec = data["srp_propagation"];
  const double t0 = ScenarioEpochTdb();
  const double Cr = sec["Cr"].get<double>();
  const double area = sec["area"].get<double>();
  const double mass = sec["mass"].get<double>();
  const double flux = sec["solar_flux"].get<double>();
  constexpr double pos_tol = 0.2;   // [m]
  constexpr double vel_tol = 1e-3;  // [m/s]

  for (const auto& kase : sec["cases"]) {
    Vec6 cart0 = ReadVec6G(kase["cart0"]);
    for (const auto& pc : kase["propagated"]) {
      double dt = pc["dt"].get<double>();
      DYNAMIC_SECTION("dt = " << dt << " s") {
        NBodyDynamics dyn;
        dyn.SetFrame(Frame::GCRF);
        dyn.AddBody(Body::Earth());
        dyn.SetUseRelativity(false);
        dyn.SetSolarFlux(flux);
        dyn.SetSrpCoefficient(Real(Cr), Real(area), Real(mass));
        Vec6 xf = Propagate(dyn, cart0, t0, dt);
        Vec6 ref = ReadVec6G(pc["cart"]);
        for (int i = 0; i < 3; i++) {
          RequireNear(xf(i), ref(i), pos_tol);
          RequireNear(xf(i + 3), ref(i + 3), vel_tol);
        }
      }
    }
  }
}
