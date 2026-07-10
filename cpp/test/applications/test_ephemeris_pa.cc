#include <lupnt/applications/ephemeris/lunanet_almanac.h>
#include <lupnt/applications/ephemeris/lunanet_ephemeris.h>
#include <lupnt/conversions/state_conversions.h>
#include <lupnt/core/constants.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "../utils.cc"

using namespace lupnt;

// These tests are self-contained (only ClassicalToCart, a pure Kepler function --
// no gravity/ephemeris data files): the MOON_PA frame handling in
// CartesianEphemeris/Almanac is an analytic rotation about +z at OMEGA_MOON, so a
// synthetic two-body arc built the same way is reproduced by the fit.
namespace {
  // Two-body arc expressed in a frame rotating about +z at `spin`: osculating
  // elements `coe_mid` (defined at the window midpoint) propagated in mean anomaly
  // and rigidly rotated by R_z(-spin*t_k). Independent construction of the model
  // that CartesianEphemeris/Almanac use for Frame::MOON_PA.
  std::pair<VecXd, MatXd> RotatingKeplerArc(const Vec6d& coe_mid, double gm, double spin, int n,
                                            double dt) {
    VecXd t_s(n);
    MatXd rv(n, 6);
    const int mid = n / 2;
    const double a = coe_mid(0);
    const double mean_motion = std::sqrt(gm / (a * a * a));
    for (int i = 0; i < n; ++i) {
      t_s(i) = dt * i;
      const double t_k = dt * (i - mid);
      Vec6 coe_k;
      coe_k << coe_mid(0), coe_mid(1), coe_mid(2), coe_mid(3), coe_mid(4),
          coe_mid(5) + mean_motion * t_k;
      const Vec6 rv_pai_r = ClassicalToCart(coe_k, Real(gm));
      const Vec6d rv_pai = rv_pai_r.cast<double>();
      const double th = -spin * t_k;
      const double c = std::cos(th), s = std::sin(th);
      Mat3d C;
      C << c, -s, 0, s, c, 0, 0, 0, 1;
      Mat3d Cdot;
      Cdot << s, c, 0, -c, s, 0, 0, 0, 0;
      Cdot *= spin;
      Vec6d rv_pa;
      rv_pa.head<3>() = C * rv_pai.head<3>();
      rv_pa.tail<3>() = C * rv_pai.tail<3>() + Cdot * rv_pai.head<3>();
      rv.row(i) = rv_pa.transpose();
    }
    return {t_s, rv};
  }
}  // namespace

TEST_CASE("applications.ephemeris.moon_pa_frame") {
  Vec6d coe;
  coe << 6.5414e6, 0.6, 56.2 * RAD, 0.2, 90.0 * RAD, 0.1;  // a, e, i, raan, argp, M
  const int n = 121;
  const double dt = 60.0;
  auto [t_s, rv_pa] = RotatingKeplerArc(coe, GM_MOON, OMEGA_MOON, n, dt);

  EphemerisFitOptions opt;
  opt.poly_order = 6;
  opt.frame = Frame::MOON_PA;
  CartesianEphemeris eph(opt);
  const VecXd params = eph.Fit(t_s, rv_pa);
  const MatXd fit = eph.Eval(t_s, params);

  // A pure rotating two-body arc is reproduced by the baseline alone (residual
  // Chebyshev coefficients ~0) in both position and velocity.
  REQUIRE((fit.leftCols(3) - rv_pa.leftCols(3)).cwiseAbs().maxCoeff() < 1e-3);
  REQUIRE((fit.rightCols(3) - rv_pa.rightCols(3)).cwiseAbs().maxCoeff() < 1e-6);

  // The osculating elements are recovered in the (inertial) PAI frame, i.e. equal
  // to the ones the arc was built from: params = [t_ref, t_fit, a, e, i, raan,
  // argp, M_ref, ...].
  REQUIRE_THAT(params(2), Catch::Matchers::WithinRel(coe(0), 1e-6));
  REQUIRE_THAT(params(3), Catch::Matchers::WithinAbs(coe(1), 1e-6));
  REQUIRE_THAT(params(4), Catch::Matchers::WithinAbs(coe(2), 1e-6));
}

TEST_CASE("applications.ephemeris.moon_ci_unchanged") {
  // spin = 0: an inertial two-body arc; the default MOON_CI frame reproduces it
  // (confirms the MOON_PA support leaves the inertial path untouched).
  Vec6d coe;
  coe << 6.5414e6, 0.6, 56.2 * RAD, 0.2, 90.0 * RAD, 0.1;
  auto [t_s, rv_ci] = RotatingKeplerArc(coe, GM_MOON, 0.0, 121, 60.0);

  EphemerisFitOptions opt;  // frame defaults to MOON_CI
  opt.poly_order = 6;
  CartesianEphemeris eph(opt);
  const MatXd fit = eph.Eval(t_s, eph.Fit(t_s, rv_ci));
  REQUIRE((fit.leftCols(3) - rv_ci.leftCols(3)).cwiseAbs().maxCoeff() < 1e-3);
  REQUIRE((fit.rightCols(3) - rv_ci.rightCols(3)).cwiseAbs().maxCoeff() < 1e-6);
}

TEST_CASE("applications.ephemeris.fourier_reduces_residual") {
  // Add a periodic once-per-orbit position wobble the Chebyshev-only model at a
  // modest order cannot absorb but a Fourier term (harmonic of the argument of
  // latitude) can.
  Vec6d coe;
  coe << 6.5414e6, 0.3, 45.0 * RAD, 0.0, 0.0, 0.0;
  const int n = 241;
  const double dt = 60.0;
  auto arc = RotatingKeplerArc(coe, GM_MOON, 0.0, n, dt);
  const VecXd t_s = arc.first;
  MatXd rv = arc.second;
  // Superimpose a small once-per-orbit oscillation in x driven by the true
  // argument of latitude (approximated by mean motion here for the perturbation).
  const double nmean = std::sqrt(GM_MOON / std::pow(coe(0), 3));
  const int mid = n / 2;
  for (int i = 0; i < n; ++i) {
    const double u = nmean * dt * (i - mid);
    rv(i, 0) += 500.0 * std::cos(u);
  }

  auto max_pos_err = [&](int num_fourier) {
    EphemerisFitOptions opt;
    opt.poly_order = 4;
    opt.num_fourier_terms = num_fourier;
    CartesianEphemeris eph(opt);
    const MatXd fit = eph.Eval(t_s, eph.Fit(t_s, rv));
    return (fit.leftCols(3) - rv.leftCols(3)).cwiseAbs().maxCoeff();
  };

  const double err_cheb = max_pos_err(0);
  const double err_fourier = max_pos_err(2);
  REQUIRE(err_fourier < err_cheb);
  REQUIRE(err_fourier < 1.0);  // Fourier captures the once-per-orbit wobble
}

TEST_CASE("applications.almanac.moon_pa_frame") {
  Vec6d coe;
  coe << 6.5414e6, 0.6, 56.2 * RAD, 0.2, 90.0 * RAD, 0.1;
  const int n = 361;
  const double dt = 600.0;
  auto [t_s, rv_pa] = RotatingKeplerArc(coe, GM_MOON, OMEGA_MOON, n, dt);

  AlmanacFitOptions opt;
  opt.poly_order = 2;
  opt.frame = Frame::MOON_PA;
  Almanac alm(opt);
  const MatXd fit = alm.Eval(t_s, alm.Fit(t_s, rv_pa));

  // A pure two-body arc has (near-)constant osculating elements plus a linear node
  // drift in the rotating frame, all captured by the low-order element polynomials.
  REQUIRE((fit.leftCols(3) - rv_pa.leftCols(3)).cwiseAbs().maxCoeff() < 10.0);
}
