#include <lupnt/conversions/anomaly_conversions.h>
#include <lupnt/conversions/state_conversions.h>
#include <lupnt/dynamics/analytical_orbit_dynamics.h>
#include <lupnt/dynamics/numerical_orbit_dynamics.h>
#include <lupnt/numerics/math_utils.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// ============================================================================
// Cross-validation against GMAT
//
// The reference values below were generated once via GMAT (R2022a) and
// checked into `cpp/test/gmat/data/gmat_reference.json` by
// `cpp/test/gmat/gen_gmat_reference.py`. This test does **not** require
// GMAT to run -- see `cpp/test/gmat/README.md` for how the fixture was
// produced and how to regenerate it.
// ============================================================================

namespace {
  Vec6 ReadVec6(const nlohmann::json& arr) {
    Vec6 v;
    for (int i = 0; i < 6; i++) v(i) = arr.at(i).get<double>();
    return v;
  }

  // Compare two angles modulo 2*pi. Some of the GMAT test cases here are
  // deliberately propagated by exactly half/1.5 orbital periods, landing
  // M/E/nu exactly on the +pi == -pi branch cut. WrapToPi(a) and WrapToPi(b)
  // can then land on opposite ends of (-pi, pi] even though `a` and `b`
  // represent the same angle, so wrap the *difference* instead (which is
  // near zero either way).
  void RequireAngleNear(Real a, Real b, double abs_error) {
    RequireNear(WrapToPi(a - b), Real(0.0), abs_error);
  }
}  // namespace

// ----------------------------------------------------------------------------
// Two-body Keplerian propagation, COE <-> Cartesian, and anomaly conversions
//
// `kepler` was generated using GMAT's point-mass-only force models
// ("TwoBodyFM" for the Earth cases, "LunaTwoBodyFM" for the lunar cases),
// so it is directly comparable to LuPNT's analytical KeplerianDynamics.
// Unlike Orekit's KeplerianPropagator (analytical, agrees with LuPNT to
// ~1e-8 m / 1e-12 rad), GMAT *numerically integrates* the point-mass
// problem (RungeKutta89, Accuracy = 1e-13), so the propagated states carry
// integrator drift: observed up to ~2e-2 m in position and ~5e-8 rad in
// the fast angles (M/E/nu) after 1.5 orbits, plus ~2e-9 rad of slow drift
// in i for the lunar cases. The propagated-state tolerances below are set
// ~5-20x above those observed integrator-drift levels (and remain many
// orders of magnitude below any real modeling error).
// ----------------------------------------------------------------------------
TEST_CASE("dynamics.kepler_gmat_reference") {
  nlohmann::json data = LoadTestJson("gmat/data/gmat_reference.json");

  for (const auto& kc : data["kepler"]) {
    Real GM = kc["gm"].get<double>();
    ClassicalOE coe0(ReadVec6(kc["coe0"]), Frame::GCRF);

    DYNAMIC_SECTION("body = " << kc["body"].get<std::string>() << ", a = " << coe0.a().val()
                              << " m, e = " << coe0.e().val()) {
      // ClassicalOE -> Cartesian
      Cart6 cart0 = ClassicalToCart(coe0, GM);
      Vec6 cart0_ref = ReadVec6(kc["cart0"]);
      for (int i = 0; i < 3; i++) {
        RequireNear(cart0(i), cart0_ref(i), 1e-3);          // [m]
        RequireNear(cart0(i + 3), cart0_ref(i + 3), 1e-6);  // [m/s]
      }

      // Cartesian -> ClassicalOE round trip. The angular elements are
      // compared modulo 2*pi (CartToClassical may return e.g. -pi/2 for an
      // input w of 3*pi/2).
      ClassicalOE coe0_rt = CartToClassical(cart0, GM);
      RequireNear(coe0_rt.a(), coe0.a(), 1e-3);
      RequireNear(coe0_rt.e(), coe0.e(), 1e-12);
      RequireNear(coe0_rt.i(), coe0.i(), 1e-12);
      RequireAngleNear(coe0_rt.Omega(), coe0.Omega(), 1e-11);
      RequireAngleNear(coe0_rt.w(), coe0.w(), 1e-11);
      RequireAngleNear(coe0_rt.M(), coe0.M(), 1e-11);

      // Anomaly conversions at epoch
      Real E0 = MeanToEccAnomaly(coe0.M(), coe0.e());
      Real nu0 = MeanToTrueAnomaly(coe0.M(), coe0.e());
      RequireAngleNear(E0, Real(kc["ecc_anomaly0"].get<double>()), 1e-9);
      RequireAngleNear(nu0, Real(kc["true_anomaly0"].get<double>()), 1e-9);

      // Orbital period
      Real period = GetOrbitalPeriod(coe0.a(), GM);
      RequireNear(period, Real(kc["period"].get<double>()), 1e-6);

      // Two-body Keplerian propagation (analytical, exact)
      KeplerianDynamics<ClassicalOE> dyn(GM);
      for (const auto& pc : kc["propagated"]) {
        Real dt = pc["dt"].get<double>();

        DYNAMIC_SECTION("dt = " << dt.val() << " s") {
          State xf = dyn.Propagate(coe0, Real(0.0), dt);
          ClassicalOE coe_f(xf);

          // Angular elements are compared modulo 2*pi (see RequireAngleNear),
          // since LuPNT and GMAT may report angles near +/-pi using
          // different (but equivalent) wrapping conventions -- some of these
          // test cases land exactly on the +pi == -pi branch cut (e.g.
          // propagation by exactly half an orbital period). The tolerances
          // absorb GMAT's RK89 integrator drift (see the TEST_CASE comment).
          Vec6 coe_ref = ReadVec6(pc["coe"]);
          RequireNear(coe_f.a(), coe_ref(0), 1e-3);
          RequireNear(coe_f.e(), coe_ref(1), 1e-9);
          RequireNear(coe_f.i(), coe_ref(2), 1e-8);
          RequireAngleNear(coe_f.Omega(), Real(coe_ref(3)), 1e-9);
          RequireAngleNear(coe_f.w(), Real(coe_ref(4)), 1e-9);
          RequireAngleNear(coe_f.M(), Real(coe_ref(5)), 1e-6);

          Cart6 cart_f = ClassicalToCart(coe_f, GM);
          Vec6 cart_ref = ReadVec6(pc["cart"]);
          for (int i = 0; i < 3; i++) {
            RequireNear(cart_f(i), cart_ref(i), 0.1);           // [m]
            RequireNear(cart_f(i + 3), cart_ref(i + 3), 1e-4);  // [m/s]
          }

          Real Ef = MeanToEccAnomaly(coe_f.M(), coe_f.e());
          Real nuf = MeanToTrueAnomaly(coe_f.M(), coe_f.e());
          RequireAngleNear(Ef, Real(pc["ecc_anomaly"].get<double>()), 1e-6);
          RequireAngleNear(nuf, Real(pc["true_anomaly"].get<double>()), 1e-6);
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------
// J2-perturbed propagation (numerical, GMAT-specific -- no Orekit equivalent)
//
// `j2_propagation` was generated using GMAT's "J2Prop"/"J2FM", a force model
// containing only the Earth point-mass + J2 term (degree/order 2,0), with the
// exact GM/R_earth/J2 values from cpp/lupnt/core/constants.h baked into a
// custom potential file (lupnt_j2_earth.cof) so both sides use identical
// constants.
//
// LuPNT's JToCartTwoBodyDynamics now evaluates the J2 acceleration in the
// body-fixed frame by default, matching GMAT's GravityField modeling
// convention for the degree-2 zonal term. The remaining tolerance primarily
// covers differences in Earth-orientation details and numerical integrators.
// ----------------------------------------------------------------------------
TEST_CASE("dynamics.j2_propagation_gmat_reference") {
  nlohmann::json data = LoadTestJson("gmat/data/gmat_reference.json");

  constexpr double pos_tol = 300.0;  // [m]
  constexpr double vel_tol = 0.15;   // [m/s]

  for (const auto& jc : data["j2_propagation"]) {
    Real GM = jc["gm"].get<double>();
    Real J2 = jc["j2"].get<double>();
    Real Re = jc["r_earth"].get<double>();

    Vec6 cart0_v = ReadVec6(jc["cart0"]);
    Cart6 cart0(cart0_v, Frame::GCRF);

    DYNAMIC_SECTION("cart0 = " << cart0_v.head<3>().transpose()) {
      JToCartTwoBodyDynamics dyn(GM, J2, Re);
      dyn.SetTimeStep(1.0);

      for (const auto& pc : jc["propagated"]) {
        Real dt = pc["dt"].get<double>();

        DYNAMIC_SECTION("dt = " << dt.val() << " s") {
          State xf = dyn.Propagate(cart0, Real(0.0), dt);
          Cart6 cart_f(xf);

          Vec6 cart_ref = ReadVec6(pc["cart"]);
          for (int i = 0; i < 3; i++) {
            RequireNear(cart_f(i), cart_ref(i), pos_tol);          // [m]
            RequireNear(cart_f(i + 3), cart_ref(i + 3), vel_tol);  // [m/s]
          }
        }
      }
    }
  }
}
