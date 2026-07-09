#include <lupnt/numerics/interpolation.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace lupnt;
using namespace Catch::Matchers;

namespace {

  TEST_CASE("numerics.interpolation.linear_1d") {
    VecXd x(4), z(4);
    x << 0, 1, 2, 3;
    z << 0, 10, 20, 30;  // z = 10 x

    REQUIRE_THAT(LinearInterp1d(x, z, 1.5), WithinAbs(15.0, 1e-9));
    REQUIRE_THAT(LinearInterp1d(x, z, 0.0), WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(LinearInterp1d(x, z, 3.0), WithinAbs(30.0, 1e-9));
    REQUIRE_THAT(LinearInterp1d(x, z, 2.25), WithinAbs(22.5, 1e-9));

    // Out of range and size mismatch throw
    REQUIRE_THROWS(LinearInterp1d(x, z, -1.0));
    REQUIRE_THROWS(LinearInterp1d(x, z, 4.0));
    VecXd z_bad(3);
    z_bad << 0, 1, 2;
    REQUIRE_THROWS(LinearInterp1d(x, z_bad, 1.0));
  }

  TEST_CASE("numerics.interpolation.linear_2d") {
    VecXd x(2), y(2);
    x << 0, 1;
    y << 0, 1;
    MatXd z(2, 2);
    z << 0, 1, 2, 3;  // z(0,0)=0, z(0,1)=1, z(1,0)=2, z(1,1)=3

    // Center is the average of the four corners
    REQUIRE_THAT(LinearInterp2d(x, y, z, 0.5, 0.5), WithinAbs(1.5, 1e-9));
    // Corners reproduce exactly
    REQUIRE_THAT(LinearInterp2d(x, y, z, 0.0, 0.0), WithinAbs(0.0, 1e-9));
    REQUIRE_THAT(LinearInterp2d(x, y, z, 1.0, 1.0), WithinAbs(3.0, 1e-9));
    // Edge midpoint along y at x=0: between z(0,0)=0 and z(0,1)=1 -> 0.5
    REQUIRE_THAT(LinearInterp2d(x, y, z, 0.0, 0.5), WithinAbs(0.5, 1e-9));

    REQUIRE_THROWS(LinearInterp2d(x, y, z, 2.0, 0.5));
    REQUIRE_THROWS(LinearInterp2d(x, y, z, 0.5, -2.0));
  }

  TEST_CASE("numerics.interpolation.lagrange") {
    VecXd x(5);
    x << 0, 1, 2, 3, 4;

    // Linear data: exact for any order
    VecXd z_lin(5);
    z_lin << 0, 2, 4, 6, 8;  // z = 2 x
    LagrangeInterpolator lin(x, 2.5, 2);
    REQUIRE_THAT(lin.Interpolate(z_lin), WithinAbs(5.0, 1e-9));

    // Quadratic data: exact with a 3-point (order 3) interpolant
    VecXd z_quad(5);
    for (int i = 0; i < 5; i++) z_quad(i) = x(i) * x(i);  // z = x^2
    LagrangeInterpolator quad(x, 2.5, 3);
    REQUIRE_THAT(quad.Interpolate(z_quad), WithinAbs(6.25, 1e-9));

    // Reproduces node values exactly
    LagrangeInterpolator at_node(x, 2.0, 3);
    REQUIRE_THAT(at_node.Interpolate(z_quad), WithinAbs(4.0, 1e-9));

    // Invalid construction / evaluation
    REQUIRE_THROWS(LagrangeInterpolator(x, 2.5, 5));   // order >= size
    REQUIRE_THROWS(LagrangeInterpolator(x, 10.0, 3));  // out of range
    VecXd z_bad(3);
    z_bad << 1, 2, 3;
    REQUIRE_THROWS(quad.Interpolate(z_bad));  // wrong size
  }

}  // namespace
