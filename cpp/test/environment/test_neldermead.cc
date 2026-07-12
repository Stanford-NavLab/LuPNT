#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "lupnt/environment/plasma/tec/neldermead.h"

using namespace pecsim;
using namespace Catch::Matchers;

TEST_CASE("environment.plasma.neldermead.paraboloid") {
  // Minimize f(x, y) = (x - 1)^2 + (y - 2)^2, min at (1, 2) with value 0.
  auto f = [](Vec2d p) { return (p[0] - 1.0) * (p[0] - 1.0) + (p[1] - 2.0) * (p[1] - 2.0); };

  NelderMead nm(f, Vec2d{-3.0, 5.0}, 1.0);
  Vec2d x = nm.minimize(1000, 1.0e-12);

  REQUIRE_THAT(x[0], WithinAbs(1.0, 1.0e-4));
  REQUIRE_THAT(x[1], WithinAbs(2.0, 1.0e-4));
  REQUIRE_THAT(nm.get_min_val(), WithinAbs(0.0, 1.0e-8));
  REQUIRE(nm.get_iter() > 0);
}

TEST_CASE("environment.plasma.neldermead.rosenbrock") {
  // Rosenbrock: f(x, y) = (1 - x)^2 + 100 (y - x^2)^2, min at (1, 1) with value 0.
  auto f = [](Vec2d p) {
    double a = 1.0 - p[0];
    double b = p[1] - p[0] * p[0];
    return a * a + 100.0 * b * b;
  };

  NelderMead nm(f, Vec2d{-1.2, 1.0}, 0.5);
  Vec2d x = nm.minimize(5000, 1.0e-14);

  REQUIRE_THAT(x[0], WithinAbs(1.0, 1.0e-2));
  REQUIRE_THAT(x[1], WithinAbs(1.0, 1.0e-2));
  REQUIRE(nm.get_min_val() < 1.0e-4);
}

TEST_CASE("environment.plasma.neldermead.fval_tol_early_stop") {
  // With a generous function-value tolerance the search returns as soon as the best vertex is
  // good enough, without exhausting the iteration budget.
  auto f = [](Vec2d p) { return p[0] * p[0] + p[1] * p[1]; };
  NelderMead nm(f, Vec2d{2.0, -2.0}, 1.0);
  Vec2d x = nm.minimize(1000, 1.0e-12, 1.0e-2);
  REQUIRE(f(x) <= 1.0e-2);
}
