#include <lupnt/numerics/math_utils.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <random>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {

  constexpr double kEps = 1e-9;

  TEST_CASE("numerics.math_utils.arange") {
    VectorX<double> a = Arange(0.0, 5.0, 1.0);
    REQUIRE(a.size() == 5);
    REQUIRE_THAT(a(0), WithinAbs(0.0, kEps));
    REQUIRE_THAT(a(4), WithinAbs(4.0, kEps));

    VectorX<int> b = Arange(2, 10, 3);  // 2, 5, 8
    REQUIRE(b.size() == 3);
    REQUIRE(b(0) == 2);
    REQUIRE(b(2) == 8);

    // Empty range (start >= stop)
    VectorX<double> c = Arange(5.0, 5.0, 1.0);
    REQUIRE(c.size() == 0);
  }

  TEST_CASE("numerics.math_utils.subsample") {
    VecX v(6);
    v << 0, 1, 2, 3, 4, 5;
    VecXd sv = Subsample(v, 2);
    REQUIRE(sv.size() == 3);
    REQUIRE_THAT(sv(0), WithinAbs(0.0, kEps));
    REQUIRE_THAT(sv(1), WithinAbs(2.0, kEps));
    REQUIRE_THAT(sv(2), WithinAbs(4.0, kEps));

    MatX m(4, 4);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) m(i, j) = 10 * i + j;
    MatXd sm = Subsample(m, 2, 2);
    REQUIRE(sm.rows() == 2);
    REQUIRE(sm.cols() == 2);
    REQUIRE_THAT(sm(1, 1), WithinAbs(22.0, kEps));  // m(2,2)
  }

  TEST_CASE("numerics.math_utils.angle_between_vecs") {
    VecX x(3), y(3);
    x << 1, 0, 0;
    y << 0, 1, 0;
    REQUIRE_THAT(AngleBetweenVecs(x, y).val(), WithinAbs(PI_OVER_TWO, 1e-12));

    VecX z(3);
    z << 1, 0, 0;
    REQUIRE_THAT(AngleBetweenVecs(x, z).val(), WithinAbs(0.0, 1e-12));

    VecX w(3);
    w << -1, 0, 0;
    REQUIRE_THAT(AngleBetweenVecs(x, w).val(), WithinAbs(PI, 1e-9));

    // Matrix (row-wise) overload
    MatX A(2, 3), B(2, 3);
    A << 1, 0, 0, 0, 1, 0;
    B << 0, 1, 0, 0, 1, 0;
    VecX ang = AngleBetweenVecs(A, B);
    REQUIRE_THAT(ang(0).val(), WithinAbs(PI_OVER_TWO, 1e-9));
    REQUIRE_THAT(ang(1).val(), WithinAbs(0.0, 1e-9));
  }

  TEST_CASE("numerics.math_utils.deg_min_sec") {
    // 30 deg 30 min 36 sec = 30.51 deg
    REQUIRE_THAT(DegMinSec2DeciDeg(30, 30, 36), WithinAbs(30.51, 1e-9));
    REQUIRE_THAT(DegMinSecToDeg(Vec3(30, 30, 36)).val(), WithinAbs(30.51, 1e-9));

    Vec3 dms = DegToDegMinSec(Real(30.51));
    REQUIRE_THAT(dms(0).val(), WithinAbs(30.0, 1e-9));
    REQUIRE_THAT(dms(1).val(), WithinAbs(30.0, 1e-9));
    REQUIRE_THAT(dms(2).val(), WithinAbs(36.0, 1e-6));
  }

  TEST_CASE("numerics.math_utils.wrap_angles") {
    // Stay off the +/-pi branch cut so the sign is unambiguous.
    REQUIRE_THAT(WrapToPi(Real(2.0 * PI + 0.5)).val(), WithinAbs(0.5, 1e-9));
    REQUIRE_THAT(WrapToPi(Real(-0.5)).val(), WithinAbs(-0.5, 1e-9));
    REQUIRE_THAT(WrapToTwoPi(Real(-0.5)).val(), WithinAbs(TWO_PI - 0.5, 1e-9));
    REQUIRE_THAT(WrapToTwoPi(Real(TWO_PI + 1.0)).val(), WithinAbs(1.0, 1e-9));

    VecX angs(2);
    angs << 2.0 * PI + 0.5, -0.5;
    VecX wp = WrapToPi(angs);
    REQUIRE_THAT(wp(0).val(), WithinAbs(0.5, 1e-9));
    REQUIRE_THAT(wp(1).val(), WithinAbs(-0.5, 1e-9));
    VecX wtp = WrapToTwoPi(angs);
    REQUIRE_THAT(wtp(1).val(), WithinAbs(TWO_PI - 0.5, 1e-9));
  }

  TEST_CASE("numerics.math_utils.decibel_roundtrip") {
    REQUIRE_THAT(DecimalToDecibel(Real(100.0)).val(), WithinAbs(20.0, 1e-9));
    REQUIRE_THAT(DecibelToDecimal(Real(20.0)).val(), WithinAbs(100.0, 1e-6));
    // Round trip
    Real x = 42.0;
    REQUIRE_THAT(DecibelToDecimal(DecimalToDecibel(x)).val(), WithinAbs(42.0, 1e-6));

    ArrX xs(2, 1);
    xs << 10.0, 1000.0;
    ArrX db = DecimalToDecibel(xs);
    REQUIRE_THAT(db(0).val(), WithinAbs(10.0, 1e-9));
    REQUIRE_THAT(db(1).val(), WithinAbs(30.0, 1e-9));
    ArrX lin = DecibelToDecimal(db);
    REQUIRE_THAT(lin(0).val(), WithinAbs(10.0, 1e-6));
    REQUIRE_THAT(lin(1).val(), WithinAbs(1000.0, 1e-6));
  }

  TEST_CASE("numerics.math_utils.min_max") {
    REQUIRE_THAT(Max(Real(3.0), Real(5.0)).val(), WithinAbs(5.0, kEps));
    REQUIRE_THAT(Min(Real(3.0), Real(5.0)).val(), WithinAbs(3.0, kEps));
    REQUIRE_THAT(MaxD(3.0, 5.0), WithinAbs(5.0, kEps));
    REQUIRE_THAT(MinD(3.0, 5.0), WithinAbs(3.0, kEps));
  }

  TEST_CASE("numerics.math_utils.rounding") {
    REQUIRE_THAT(round(Real(2.567), 1).val(), WithinAbs(2.6, 1e-9));
    REQUIRE_THAT(round(Real(2.4)).val(), WithinAbs(2.0, 1e-9));
    REQUIRE_THAT(frac(Real(2.75)).val(), WithinAbs(0.75, 1e-9));
    REQUIRE_THAT(ceil(Real(2.1)).val(), WithinAbs(3.0, 1e-9));
    REQUIRE_THAT(floor(Real(2.9)).val(), WithinAbs(2.0, 1e-9));
    REQUIRE_THAT(mod(Real(7.0), Real(3.0)).val(), WithinAbs(1.0, 1e-9));
  }

  TEST_CASE("numerics.math_utils.trig_degrees") {
    REQUIRE_THAT(sind(Real(30.0)).val(), WithinAbs(0.5, 1e-12));
    REQUIRE_THAT(cosd(Real(60.0)).val(), WithinAbs(0.5, 1e-12));
    REQUIRE_THAT(tand(Real(45.0)).val(), WithinAbs(1.0, 1e-12));
  }

  TEST_CASE("numerics.math_utils.safe_trig") {
    // At the exact boundary (and just past it, as from round-off) the result stays
    // finite instead of producing NaN.
    REQUIRE_THAT(safe_acos(Real(1.0)).val(), WithinAbs(0.0, 1e-6));
    REQUIRE_THAT(safe_acos(Real(-1.0)).val(), WithinAbs(PI, 1e-6));
    REQUIRE_THAT(safe_acos(Real(0.0)).val(), WithinAbs(PI_OVER_TWO, 1e-9));
    REQUIRE_THAT(safe_asin(Real(1.0)).val(), WithinAbs(PI_OVER_TWO, 1e-6));
    REQUIRE_THAT(safe_asin(Real(-1.0)).val(), WithinAbs(-PI_OVER_TWO, 1e-6));
    REQUIRE_THAT(safe_asin(Real(0.0)).val(), WithinAbs(0.0, 1e-9));
    // Ordinary in-range values match the standard functions.
    REQUIRE_THAT(safe_acos(Real(0.5)).val(), WithinAbs(std::acos(0.5), 1e-12));
    REQUIRE_THAT(safe_asin(Real(0.5)).val(), WithinAbs(std::asin(0.5), 1e-12));
  }

  TEST_CASE("numerics.math_utils.statistics") {
    VecXd x(4);
    x << 1.0, 2.0, 3.0, 4.0;
    REQUIRE_THAT(RootMeanSquare(x), WithinAbs(std::sqrt(30.0 / 4.0), 1e-9));
    // sample std of {1,2,3,4} = sqrt(5/3)
    REQUIRE_THAT(Std(x), WithinAbs(std::sqrt(5.0 / 3.0), 1e-9));
    // percentile index = ceil(p*(n-1))
    REQUIRE_THAT(Percentile(x, 0.0), WithinAbs(1.0, 1e-9));
    REQUIRE_THAT(Percentile(x, 1.0), WithinAbs(4.0, 1e-9));
  }

  TEST_CASE("numerics.math_utils.rotations") {
    // Passive rotation of pi/2 about z takes x-hat to -y-hat
    Vec3 xhat(1, 0, 0);
    Vec3 r = RotZ(Real(PI_OVER_TWO)) * xhat;
    REQUIRE_THAT(r(0).val(), WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(r(1).val(), WithinAbs(-1.0, 1e-12));

    // Rotation matrices are orthonormal (R^T R = I)
    Mat3 Rx = RotX(Real(0.3));
    Mat3 Ry = RotY(Real(-0.7));
    Mat3 Rz = RotZ(Real(1.1));
    RequireNear(Rx.transpose() * Rx, Mat3::Identity(), 1e-12);
    RequireNear(Ry.transpose() * Ry, Mat3::Identity(), 1e-12);
    RequireNear(Rz.transpose() * Rz, Mat3::Identity(), 1e-12);

    // RotationAngle recovers the rotation magnitude
    REQUIRE_THAT(RotationAngle(RotZ(Real(1.1))).val(), WithinAbs(1.1, 1e-9));
  }

  TEST_CASE("numerics.math_utils.rot_dot") {
    // d/dt RotZ(theta(t)) numerically matches RotZdot for a given angle rate
    Real theta = 0.4;
    Real dtheta = 1.3;
    double h = 1e-6;
    Mat3 num = (RotZ(theta + dtheta * h) - RotZ(theta - dtheta * h)) / (2.0 * h);
    RequireNear(RotZdot(theta, dtheta), num, 1e-5);

    Mat3 numx = (RotX(theta + dtheta * h) - RotX(theta - dtheta * h)) / (2.0 * h);
    RequireNear(RotXdot(theta, dtheta), numx, 1e-5);
    Mat3 numy = (RotY(theta + dtheta * h) - RotY(theta - dtheta * h)) / (2.0 * h);
    RequireNear(RotYdot(theta, dtheta), numy, 1e-5);
  }

  TEST_CASE("numerics.math_utils.skew") {
    Vec3 a(1, 2, 3);
    Vec3 b(4, 5, 6);
    Vec3 sk = Skew(a) * b;
    Vec3 cr = a.cross(b);
    RequireNear(sk, cr, 1e-12);
    // Skew matrix is antisymmetric
    Mat3 S = Skew(a);
    RequireNear(S + S.transpose(), Mat3::Zero(), 1e-12);
  }

  TEST_CASE("numerics.math_utils.eigen_to_std") {
    VecXd v(3);
    v << 1.5, 2.5, 3.5;
    std::vector<double> sv = EigenToVector<double>(v);
    REQUIRE(sv.size() == 3);
    REQUIRE_THAT(sv[1], WithinAbs(2.5, 1e-9));
  }

  TEST_CASE("numerics.math_utils.block_diagonal") {
    Matrix<Real, 2, 2> A;
    A << 1, 2, 3, 4;
    Matrix<Real, 2, 2> B;
    B << 5, 6, 7, 8;
    Matrix<Real, 4, 4> C = BlockDiagonal(A, B);
    REQUIRE_THAT(C(0, 0).val(), WithinAbs(1.0, kEps));
    REQUIRE_THAT(C(3, 3).val(), WithinAbs(8.0, kEps));
    // Off-diagonal blocks are zero
    REQUIRE_THAT(C(0, 2).val(), WithinAbs(0.0, kEps));
    REQUIRE_THAT(C(2, 0).val(), WithinAbs(0.0, kEps));
  }

  TEST_CASE("numerics.math_utils.vech") {
    MatX m(3, 3);
    m << 1, 99, 99, 2, 3, 99, 4, 5, 6;  // only lower triangle matters
    VecX v = Vech(m);
    REQUIRE(v.size() == 6);
    double expected[6] = {1, 2, 3, 4, 5, 6};
    for (int i = 0; i < 6; i++) REQUIRE_THAT(v(i).val(), WithinAbs(expected[i], kEps));
  }

  TEST_CASE("numerics.math_utils.unpack") {
    Vec3 v(7.0, 8.0, 9.0);
    auto [a, b, c] = Unpack(Vec3(v));
    REQUIRE_THAT(a.val(), WithinAbs(7.0, kEps));
    REQUIRE_THAT(b.val(), WithinAbs(8.0, kEps));
    REQUIRE_THAT(c.val(), WithinAbs(9.0, kEps));
  }

  TEST_CASE("numerics.math_utils.linear_solvers") {
    MatXd A(2, 2);
    A << 2, 0, 0, 4;
    VecXd b(2);
    b << 6, 8;
    VecXd x = SolveLinearEqSVD(A, b);
    REQUIRE_THAT(x(0), WithinAbs(3.0, 1e-9));
    REQUIRE_THAT(x(1), WithinAbs(2.0, 1e-9));

    // Pseudo-inverse of an invertible matrix is its inverse
    MatXd Ainv = PseudoInverse(A);
    MatXd I = A * Ainv;
    REQUIRE_THAT(I(0, 0), WithinAbs(1.0, 1e-9));
    REQUIRE_THAT(I(1, 1), WithinAbs(1.0, 1e-9));
    REQUIRE_THAT(I(0, 1), WithinAbs(0.0, 1e-9));

    // Matrix RHS overload
    MatXd Bm(2, 2);
    Bm << 6, 2, 8, 4;
    MatXd X = SolveLinearEqSVD(A, Bm);
    REQUIRE_THAT(X(0, 0), WithinAbs(3.0, 1e-9));
    REQUIRE_THAT(X(1, 1), WithinAbs(1.0, 1e-9));
  }

  TEST_CASE("numerics.math_utils.jacobian_parallel") {
    // f(x) = [2*x0 + x1, x0 - 3*x1]; Jacobian is constant [[2,1],[1,-3]]
    std::function<VecX(const VecX&)> f = [](const VecX& x) {
      VecX y(2);
      y(0) = 2.0 * x(0) + x(1);
      y(1) = x(0) - 3.0 * x(1);
      return y;
    };
    VecX x0(2);
    x0 << 1.0, 2.0;
    MatXd J;
    VecX y = JacobianParallel(f, x0, J);
    REQUIRE_THAT(y(0).val(), WithinAbs(4.0, 1e-9));
    REQUIRE_THAT(y(1).val(), WithinAbs(-5.0, 1e-9));
    REQUIRE(J.rows() == 2);
    REQUIRE(J.cols() == 2);
    REQUIRE_THAT(J(0, 0), WithinAbs(2.0, 1e-6));
    REQUIRE_THAT(J(0, 1), WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(J(1, 0), WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(J(1, 1), WithinAbs(-3.0, 1e-6));
  }

  TEST_CASE("numerics.math_utils.sample_normal") {
    std::mt19937 rng(12345);
    // Empirical mean of many scalar samples is close to the requested mean
    int n = 20000;
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += SampleNormal(1.0, 2.0, &rng).val();
    REQUIRE_THAT(sum / n, WithinAbs(1.0, 0.1));

    // Multivariate: identity covariance, empirical mean near requested mean
    VecX mean(2);
    mean << 3.0, -1.0;
    MatX cov = MatX::Identity(2, 2);
    std::mt19937 rng2(999);
    MatX samples = SampleMvNormal(mean, cov, 20000, &rng2);
    REQUIRE(samples.rows() == 20000);
    REQUIRE(samples.cols() == 2);
    REQUIRE_THAT(samples.col(0).cast<double>().mean(), WithinAbs(3.0, 0.1));
    REQUIRE_THAT(samples.col(1).cast<double>().mean(), WithinAbs(-1.0, 0.1));
  }

  // Reference values from standard tables (Abramowitz & Stegun / SciPy).
  TEST_CASE("numerics.math_utils.special_functions") {
    SECTION("erfc matches known complementary-error-function values") {
      REQUIRE_THAT(erfc(Real(0.0)).val(), WithinAbs(1.0, 1e-12));
      REQUIRE_THAT(erfc(Real(1.0)).val(), WithinAbs(0.15729920705028513, 1e-9));
      REQUIRE_THAT(erfc(Real(2.0)).val(), WithinAbs(0.004677734981047266, 1e-9));
      // erfc(-x) = 2 - erfc(x).
      REQUIRE_THAT(erfc(Real(-1.0)).val(), WithinAbs(2.0 - 0.15729920705028513, 1e-9));
    }

    SECTION("qfunc is the upper-tail standard-normal probability") {
      REQUIRE_THAT(qfunc(Real(0.0)).val(), WithinAbs(0.5, 1e-12));
      REQUIRE_THAT(qfunc(Real(1.0)).val(), WithinAbs(0.15865525393145707, 1e-9));
      REQUIRE_THAT(qfunc(Real(2.0)).val(), WithinAbs(0.022750131948179195, 1e-9));
      // Q(-x) = 1 - Q(x).
      REQUIRE_THAT(qfunc(Real(-1.5)).val(), WithinAbs(1.0 - qfunc(Real(1.5)).val(), 1e-12));
    }

    SECTION("J0Bessel matches Bessel-J0 reference values") {
      REQUIRE_THAT(J0Bessel(Real(0.0)).val(), WithinAbs(1.0, 1e-12));
      REQUIRE_THAT(J0Bessel(Real(1.0)).val(), WithinAbs(0.7651976865579666, 1e-8));
      REQUIRE_THAT(J0Bessel(Real(2.0)).val(), WithinAbs(0.22389077914123567, 1e-7));
    }

    SECTION("J1Bessel matches Bessel-J1 reference values") {
      REQUIRE_THAT(J1Bessel(Real(0.0)).val(), WithinAbs(0.0, 1e-12));
      REQUIRE_THAT(J1Bessel(Real(1.0)).val(), WithinAbs(0.4400505857449335, 1e-8));
      REQUIRE_THAT(J1Bessel(Real(2.0)).val(), WithinAbs(0.5767248077568734, 1e-7));
    }
  }

}  // namespace
