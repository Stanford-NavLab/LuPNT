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
// Cross-validation against Orekit
//
// The reference values below were generated once via Orekit and checked
// into `cpp/test/orekit/data/orekit_reference.json` by
// `cpp/test/orekit/gen_orekit_reference.py`. This test does **not**
// require Orekit/Java to run -- see `cpp/test/orekit/README.md` for how
// the fixture was produced and how to regenerate it.
// ============================================================================

namespace {
  Vec6 ReadVec6(const nlohmann::json& arr) {
    Vec6 v;
    for (int i = 0; i < 6; i++) v(i) = arr.at(i).get<double>();
    return v;
  }
  Vec3 ReadVec3Json(const nlohmann::json& arr) {
    Vec3 v;
    for (int i = 0; i < 3; i++) v(i) = arr.at(i).get<double>();
    return v;
  }

  // Compare two angles modulo 2*pi. Some of the test cases here are
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
// The `kepler` cases span LEO/MEO-HEO/GEO/Molniya Earth orbits (GM_EARTH)
// plus low-lunar-orbit and ELFO cases (GM_MOON); the `gm` field of each case
// selects the gravitational parameter, so the same code exercises both
// central bodies.
// ----------------------------------------------------------------------------
TEST_CASE("dynamics.kepler_orekit_reference") {
  nlohmann::json data = LoadTestJson("orekit/data/orekit_reference.json");

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

          // Angular elements are compared modulo 2*pi (see RequireAngleNear):
          // LuPNT and Orekit may report angles near +/-pi (or 3*pi/2 vs
          // -pi/2) using different but equivalent wrapping conventions, and
          // some cases land exactly on the branch cut (e.g. propagation by
          // exactly half an orbital period from M0 = 0).
          Vec6 coe_ref = ReadVec6(pc["coe"]);
          RequireNear(coe_f.a(), coe_ref(0), 1e-3);
          RequireNear(coe_f.e(), coe_ref(1), 1e-9);
          RequireNear(coe_f.i(), coe_ref(2), 1e-12);
          RequireAngleNear(coe_f.Omega(), Real(coe_ref(3)), 1e-12);
          RequireAngleNear(coe_f.w(), Real(coe_ref(4)), 1e-9);
          RequireAngleNear(coe_f.M(), Real(coe_ref(5)), 1e-9);

          Cart6 cart_f = ClassicalToCart(coe_f, GM);
          Vec6 cart_ref = ReadVec6(pc["cart"]);
          for (int i = 0; i < 3; i++) {
            RequireNear(cart_f(i), cart_ref(i), 1e-2);          // [m]
            RequireNear(cart_f(i + 3), cart_ref(i + 3), 1e-5);  // [m/s]
          }

          Real Ef = MeanToEccAnomaly(coe_f.M(), coe_f.e());
          Real nuf = MeanToTrueAnomaly(coe_f.M(), coe_f.e());
          RequireAngleNear(Ef, Real(pc["ecc_anomaly"].get<double>()), 1e-9);
          RequireAngleNear(nuf, Real(pc["true_anomaly"].get<double>()), 1e-9);
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------
// J2 perturbing acceleration (formula-level, frame-independent check)
// ----------------------------------------------------------------------------
TEST_CASE("dynamics.j2_acceleration_orekit_reference") {
  nlohmann::json data = LoadTestJson("orekit/data/orekit_reference.json");

  JToCartTwoBodyDynamics dyn(GM_EARTH, J2_EARTH, R_EARTH, Frame::GCRF, Frame::GCRF);

  for (const auto& tc : data["j2_acceleration"]) {
    Vec3 r = ReadVec3Json(tc["r"]);
    Vec3 a_j2_ref = ReadVec3Json(tc["a_j2"]);

    DYNAMIC_SECTION("r = " << r.transpose()) {
      Cart6 x0(r, Vec3::Zero(), Frame::GCRF);
      VecX rv_dot = dyn.ComputeRates(Real(0.0), x0);

      Real r_norm = r.norm();
      Vec3 a_two_body = -GM_EARTH * r / pow(r_norm, 3);
      Vec3 a_j2 = rv_dot.tail<3>() - a_two_body;

      for (int i = 0; i < 3; i++) RequireNear(a_j2(i), a_j2_ref(i), 1e-12);
    }
  }
}

// ----------------------------------------------------------------------------
// J2-perturbed numerical propagation
//
// The Orekit reference was generated with the J2OnlyPerturbation force model
// evaluated in the **GCRF (inertial)** frame. This test configures
// LuPNT's JToCartTwoBodyDynamics with Frame::GCRF as its "body-fixed" frame
// to preserve that inertial-axis reference case. Unlike the
// GMAT j2_propagation comparison (where GMAT evaluates J2 in the body-fixed
// frame, an expected ~10-200 m modeling difference), this comparison
// isolates pure *numerical integration* differences: LuPNT's fixed-step
// RK4 (1 s) vs Orekit's adaptive DormandPrince853 (1e-6 m position
// tolerance). The observed agreement is ~1e-5 m / 1e-8 m/s over up to 1.5
// orbits across all three cases (MEO-HEO, LEO, Molniya), so the tolerances
// below are ~100x the observed differences while remaining ~6 orders of
// magnitude tighter than a force-model error (km scale) and ~5 orders
// tighter than the body-fixed-vs-inertial J2 modeling difference seen in
// the GMAT comparison (tens to hundreds of m).
// ----------------------------------------------------------------------------
TEST_CASE("dynamics.j2_propagation_orekit_reference") {
  nlohmann::json data = LoadTestJson("orekit/data/orekit_reference.json");

  constexpr double pos_tol = 1e-3;  // [m]
  constexpr double vel_tol = 1e-6;  // [m/s]

  for (const auto& jc : data["j2_propagation"]) {
    Real GM = jc["gm"].get<double>();
    Real J2 = jc["j2"].get<double>();
    Real Re = jc["r_earth"].get<double>();

    Vec6 cart0_v = ReadVec6(jc["cart0"]);
    Cart6 cart0(cart0_v, Frame::GCRF);

    DYNAMIC_SECTION("cart0 = " << cart0_v.head<3>().transpose()) {
      JToCartTwoBodyDynamics dyn(GM, J2, Re, Frame::GCRF, Frame::GCRF);
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
