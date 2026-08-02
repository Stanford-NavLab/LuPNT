#include <lupnt/core/constants.h>
#include <lupnt/environment/body.h>
#include <lupnt/environment/forces.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

// ============================================================================
// Cross-validation of individual force models against Orekit.
//
// Unlike the propagation-level `dynamics.*_orekit_reference` tests, these are
// formula-level acceleration checks (like `dynamics.j2_acceleration_orekit_...`):
// they feed LuPNT's acceleration functions the exact same inputs Orekit used
// and compare the acceleration vector directly, so there is no integrator or
// Earth-orientation noise -- only the force-model algorithm is exercised.
//
// The reference values are in `cpp/test/orekit/data/orekit_reference.json`,
// produced by `gen_orekit_reference.py` from Orekit's HolmesFeatherstone,
// ThirdBodyAttraction and SolarRadiationPressure models. See
// `cpp/test/orekit/README.md`.
// ============================================================================

namespace {
  Vec3 ReadVec3J(const nlohmann::json& arr) {
    Vec3 v;
    for (int i = 0; i < 3; i++) v(i) = arr.at(i).get<double>();
    return v;
  }
}  // namespace

// ----------------------------------------------------------------------------
// Spherical-harmonic gravity (Earth 8x8 EGM96, Moon 12x12 GRGM1200B)
//
// Both sides use the SAME `.cof` coefficient file: the generator parses it into
// an Orekit ICGEM field, and here LuPNT loads it via ReadHarmonicGravityField.
// The comparison therefore isolates the Cunningham/Pines recursion and the
// coefficient normalization convention, not the adopted coefficients. Observed
// agreement is ~5e-15 m/s^2 (machine precision); the tolerance is far above
// that yet ~10 orders below any real force-model error.
// ----------------------------------------------------------------------------
TEST_CASE("dynamics.gravity_field_orekit_reference") {
  nlohmann::json data = LoadTestJson("orekit/data/orekit_reference.json");
  constexpr double tol = 1e-9;  // [m/s^2]

  for (auto it = data["gravity_acceleration"].begin(); it != data["gravity_acceleration"].end();
       ++it) {
    const std::string body = it.key();
    const auto& sec = it.value();
    const std::string cof = sec["cof_file"].get<std::string>();
    const int n = sec["n_max"].get<int>();
    const int m = sec["m_max"].get<int>();

    DYNAMIC_SECTION("body = " << body << " (" << cof << ", " << n << "x" << m << ")") {
      GravityField<Real> gf = ReadHarmonicGravityField<Real>(cof, n, m, /*normalized=*/true);
      // The reference stored the file's GM/R; confirm LuPNT loaded the same.
      RequireNear(gf.GM, Real(sec["GM"].get<double>()), 1e-3);
      RequireNear(gf.R, Real(sec["R"].get<double>()), 1e-6);

      for (const auto& c : sec["cases"]) {
        Vec3 r_bf = ReadVec3J(c["r_bf"]);
        Vec3 a_ref = ReadVec3J(c["a_bf"]);
        Vec3 a = AccelarationGravityField<Real>(r_bf, gf.GM, gf.R, gf.CS, n, m);
        for (int i = 0; i < 3; i++) RequireNear(a(i), a_ref(i), tol);
      }
    }
  }
}

// ----------------------------------------------------------------------------
// Third-body point-mass perturbation (Sun and Moon on an Earth orbit)
//
// LuPNT's AccelerationPointMass is fed the perturber position Orekit used (its
// DE440 lookup, stored per case), so the ephemeris is decoupled from the
// perturbation formula -- the DE440 agreement itself is checked by the
// `data.ephemeris_orekit_reference` test. Observed agreement ~2e-18 m/s^2.
// ----------------------------------------------------------------------------
TEST_CASE("dynamics.third_body_orekit_reference") {
  nlohmann::json data = LoadTestJson("orekit/data/orekit_reference.json");
  constexpr double tol = 1e-12;  // [m/s^2]

  for (const auto& c : data["third_body_acceleration"]["EARTH_CENTERED"]) {
    Vec3 r = ReadVec3J(c["r"]);
    DYNAMIC_SECTION("r = " << r.transpose()) {
      for (const auto& p : c["perturbers"]) {
        Real GM = p["GM"].get<double>();
        Vec3 s = ReadVec3J(p["s"]);
        Vec3 a_ref = ReadVec3J(p["a"]);
        DYNAMIC_SECTION("perturber = " << p["body"].get<std::string>()) {
          Vec3 a = AccelerationPointMass(r, s, GM);
          for (int i = 0; i < 3; i++) RequireNear(a(i), a_ref(i), tol);
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------
// Solar radiation pressure, cannonball, fully sunlit
//
// Orekit's reference pressure at 1 AU is its own adopted constant; the fixture
// stores it (`P0`) so LuPNT uses the same value -- the check isolates the
// cannonball formula and its Sun->satellite direction / inverse-square /
// AU^2 convention, not the flux. Observed agreement ~7e-24 m/s^2.
// ----------------------------------------------------------------------------
TEST_CASE("dynamics.srp_orekit_reference") {
  nlohmann::json data = LoadTestJson("orekit/data/orekit_reference.json");
  constexpr double tol = 1e-15;  // [m/s^2]

  for (const auto& c : data["srp_acceleration"]["cases"]) {
    Vec3 r = ReadVec3J(c["r"]);
    Vec3 r_sun = ReadVec3J(c["r_sun"]);
    Real bcoeff = c["Cr"].get<double>() * c["area"].get<double>() / c["mass"].get<double>();
    Real P0 = c["P0"].get<double>();
    Real AU_ref = c["AU"].get<double>();
    Vec3 a_ref = ReadVec3J(c["a"]);

    DYNAMIC_SECTION("r = " << r.transpose()) {
      Vec3 a = AccelerationSolarRadiation(r, r_sun, bcoeff, P0, AU_ref);
      for (int i = 0; i < 3; i++) RequireNear(a(i), a_ref(i), tol);
    }
  }
}
