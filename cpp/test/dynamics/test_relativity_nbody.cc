// Validation for the n-body relativistic perturbative acceleration
// (AccelerationRelativisticNBody), which implements Moyer (2000), "Formulation
// for Observed and Computed Values of DSN Data Types for Navigation" (JPL Pub
// 00-7), Eq. (4-26) with the leading Newtonian term removed.
//
// Key check: Moyer states (p. 4-43) that Eq. (4-26), simplified to a single
// perturbing body and with the Newtonian term removed, reduces exactly to the
// 1-body Schwarzschild form Eq. (4-61):
//   a = (mu/(c^2 r^3)) { [2(beta+gamma) mu/r - gamma |v|^2] r + 2(1+gamma)(r.v) v }
// We use that closed form as an independent reference.

#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <vector>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // Moyer Eq. (4-61): 1-body Schwarzschild relativistic perturbative acceleration
  // of a spacecraft about a body at rest at the origin.
  Vec3d Eq461(const Vec3d& r, const Vec3d& v, double mu, double c, double beta, double gamma) {
    double rn = r.norm();
    double v2 = v.squaredNorm();
    double rv = r.dot(v);
    return mu / (c * c * std::pow(rn, 3))
           * ((2.0 * (beta + gamma) * mu / rn - gamma * v2) * r + 2.0 * (1.0 + gamma) * rv * v);
  }
}  // namespace

TEST_CASE("dynamics.relativity_nbody.single_body_reduces_to_schwarzschild") {
  // SI units throughout (m, m/s, m^3/s^2).
  const double mu = 3.986004418e14;  // Earth GM
  const double c = 299792458.0;
  const Vec3 r_sc(7.0e6, 1.0e6, 2.0e6);
  const Vec3 v_sc(1.0e3, 7.0e3, 0.5e3);

  // A single perturbing body at rest at the barycentric origin.
  std::vector<Vec3> r_bodies = {Vec3::Zero()};
  std::vector<Vec3> v_bodies = {Vec3::Zero()};
  std::vector<Real> mu_bodies = {Real(mu)};

  Vec3 a = AccelerationRelativisticNBody(r_sc, v_sc, r_bodies, v_bodies, mu_bodies, Real(c));
  Vec3d ref = Eq461(r_sc.cast<double>(), v_sc.cast<double>(), mu, c, 1.0, 1.0);
  for (int i = 0; i < 3; ++i) REQUIRE_THAT(a(i).val(), WithinRel(ref(i), 1e-9));

  // Sanity: the correction is tiny relative to Newtonian gravity (~10^-9 or less).
  double a_newton = mu / r_sc.squaredNorm().val();
  REQUIRE(a.norm().val() < 1e-7 * a_newton);
}

TEST_CASE("dynamics.relativity_nbody.ppn_parameters_scale_terms") {
  // With gamma = beta = 0 the velocity-squared and cross terms vanish or change
  // sign relative to GR; check the result stays finite and differs from GR.
  const double mu = 3.986004418e14, c = 299792458.0;
  const Vec3 r_sc(7.0e6, 0.0, 0.0), v_sc(0.0, 7.5e3, 0.0);
  std::vector<Vec3> r_bodies = {Vec3::Zero()};
  std::vector<Vec3> v_bodies = {Vec3::Zero()};
  std::vector<Real> mu_bodies = {Real(mu)};

  Vec3 a_gr = AccelerationRelativisticNBody(r_sc, v_sc, r_bodies, v_bodies, mu_bodies, Real(c),
                                            Real(1.0), Real(1.0));
  Vec3 a_pp = AccelerationRelativisticNBody(r_sc, v_sc, r_bodies, v_bodies, mu_bodies, Real(c),
                                            Real(0.0), Real(0.0));
  Vec3d ref_gr = Eq461(r_sc.cast<double>(), v_sc.cast<double>(), mu, c, 1.0, 1.0);
  Vec3d ref_pp = Eq461(r_sc.cast<double>(), v_sc.cast<double>(), mu, c, 0.0, 0.0);
  for (int i = 0; i < 3; ++i) {
    REQUIRE(std::isfinite(a_pp(i).val()));
    REQUIRE_THAT(a_gr(i).val(), WithinRel(ref_gr(i), 1e-9));
    REQUIRE_THAT(a_pp(i).val(), WithinRel(ref_pp(i), 1e-9));
  }
  REQUIRE((a_gr - a_pp).norm().val() > 0.0);
}

TEST_CASE("dynamics.relativity_nbody.multibody_finite_and_small") {
  // Sun + Earth + a Moon-like body: a full n-body evaluation must be finite,
  // non-NaN, and small relative to the central body's Newtonian acceleration.
  const double c = 299792458.0;
  const Vec3 r_sc(7.0e6, 0.0, 0.0);
  const Vec3 v_sc(0.0, 7.5e3, 0.0);
  std::vector<Vec3> r_bodies = {Vec3::Zero(), Vec3(1.5e11, 0.0, 0.0), Vec3(3.8e8, 0.0, 0.0)};
  std::vector<Vec3> v_bodies = {Vec3::Zero(), Vec3(0.0, 3.0e4, 0.0), Vec3(0.0, 1.0e3, 0.0)};
  std::vector<Real> mu_bodies = {Real(3.986004418e14), Real(1.32712440018e20), Real(4.9028e12)};

  Vec3 a = AccelerationRelativisticNBody(r_sc, v_sc, r_bodies, v_bodies, mu_bodies, Real(c));
  for (int i = 0; i < 3; ++i) REQUIRE(std::isfinite(a(i).val()));
  double a_newton = 3.986004418e14 / r_sc.squaredNorm().val();
  REQUIRE(a.norm().val() > 0.0);
  REQUIRE(a.norm().val() < 1e-6 * a_newton);
}
