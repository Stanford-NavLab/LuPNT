/**
 * @file nelder_mead.cpp
 * @author Keidai Iiyama
 * @brief Implementation of the Nelder-Mead optimization algorithm in 2D (Used
 * to correct intial angles in ray tracing)
 * @version 0.1
 * @date 2025-02-17
 */

#include "lupnt/environment/plasma/tec/neldermead.h"

namespace pecsim {

  NelderMead::NelderMead(std::function<double(Vec2d)> func) { f = func; }

  NelderMead::NelderMead(std::function<double(Vec2d)> func, Vec2d x0, double step, bool debug) {
    f = func;
    set_initial_simplex(x0, step, debug);
  }

  void NelderMead::set_initial_simplex(Vec2d x0, double step, bool debug) {
    // Initialize 2D simplex (3 vertices)
    Vec2d x1 = {0, 0};
    Vec2d x2 = {1, 0};
    Vec2d x3 = {0.5, std::sqrt(3) / 2.0};

    // Compute centroid
    Vec2d centroid = {(x1[0] + x2[0] + x3[0]) / 3.0, (x1[1] + x2[1] + x3[1]) / 3.0};

    // Offset points so that centroid is at x0
    for (int i = 0; i < 2; ++i) {
      x1[i] = x0[i] + step * (x1[i] - centroid[i]);
      x2[i] = x0[i] + step * (x2[i] - centroid[i]);
      x3[i] = x0[i] + step * (x3[i] - centroid[i]);
    }

    if (debug) {
      std::cout << "Initial simplex vertices:\n"
                << "  x1: (" << x1[0] << ", " << x1[1] << ")\n"
                << "  x2: (" << x2[0] << ", " << x2[1] << ")\n"
                << "  x3: (" << x3[0] << ", " << x3[1] << ")\n";
    }

    simplex[0] = {x1, f(x1)};
    simplex[1] = {x2, f(x2)};
    simplex[2] = {x3, f(x3)};
  }

  void NelderMead::set_initial_simplex(Vec2d x0, double f0, double step, bool debug) {
    // Initialize 2D simplex (3 vertices)
    Vec2d x1 = x0;                   // Start at x0
    Vec2d x2 = x0 + Vec2d{step, 0};  // Step right from x0
    Vec2d x3 = x0 + Vec2d{0, step};  // Step up from x0

    if (debug) {
      std::cout << "Initial simplex vertices:\n"
                << "  x1: (" << x1[0] << ", " << x1[1] << ")\n"
                << "  x2: (" << x2[0] << ", " << x2[1] << ")\n"
                << "  x3: (" << x3[0] << ", " << x3[1] << ")\n";
    }

    simplex[0] = {x1, f0};  // use provided function value for x0
    simplex[1] = {x2, f(x2)};
    simplex[2] = {x3, f(x3)};
  }

  Vec2d NelderMead::minimize(int max_iters, double tol, double fval_tol, bool debug) {
    if (debug) {
      std::cout << "Starting Nelder-Mead optimization with initial points:\n";
      for (const auto& vertex : simplex) {
        std::cout << "  Point: (" << vertex.x[0] << ", " << vertex.x[1] << ") fval: " << vertex.fval
                  << "\n";
      }
    }

    for (int iter = 0; iter < max_iters; ++iter) {
      iter_ = iter;
      // Sort simplex: best -> worst
      std::sort(simplex.begin(), simplex.end(),
                [](const SimplexVertex& a, const SimplexVertex& b) { return a.fval < b.fval; });

      min_val_ = simplex[0].fval;  // Update minimum value found

      if (debug) {
        std::cout << "Iteration " << iter << ": "
                  << "Best: (" << simplex[0].x[0] << ", " << simplex[0].x[1] << ") "
                  << "fval: " << simplex[0].fval << " | "
                  << "Second: (" << simplex[1].x[0] << ", " << simplex[1].x[1] << ") "
                  << "fval: " << simplex[1].fval << " | "
                  << "Worst: (" << simplex[2].x[0] << ", " << simplex[2].x[1] << ") "
                  << "fval: " << simplex[2].fval << "\n";
      }

      // Check convergence
      double spread = std::abs(simplex[2].fval - simplex[0].fval);
      if ((spread < tol) || (simplex[0].fval <= fval_tol)) return simplex[0].x;

      Vec2d centroid = mean_point({simplex[0].x, simplex[1].x});

      // Reflection
      Vec2d xr = reflect(centroid, simplex[2].x);
      double fr = f(xr);

      if (fr < simplex[0].fval) {
        // Expansion
        Vec2d xe = expand(centroid, xr);
        double fe = f(xe);
        if (fe < fr)
          simplex[2] = {xe, fe};
        else
          simplex[2] = {xr, fr};
      } else if (fr < simplex[1].fval) {
        simplex[2] = {xr, fr};
      } else {
        // Contraction
        Vec2d xc = contract(centroid, simplex[2].x);
        double fc = f(xc);
        if (fc < simplex[2].fval)
          simplex[2] = {xc, fc};
        else
          shrink();
      }
    }
    return simplex[0].x;
  }

  Vec2d NelderMead::mean_point(const std::array<Vec2d, 2>& points) {
    return {(points[0][0] + points[1][0]) / 2.0, (points[0][1] + points[1][1]) / 2.0};
  }

  Vec2d NelderMead::reflect(const Vec2d& c, const Vec2d& worst, double alpha) {
    return {c[0] + alpha * (c[0] - worst[0]), c[1] + alpha * (c[1] - worst[1])};
  }

  Vec2d NelderMead::expand(const Vec2d& c, const Vec2d& reflected, double gamma) {
    return {c[0] + gamma * (reflected[0] - c[0]), c[1] + gamma * (reflected[1] - c[1])};
  }

  Vec2d NelderMead::contract(const Vec2d& c, const Vec2d& worst, double rho) {
    return {c[0] + rho * (worst[0] - c[0]), c[1] + rho * (worst[1] - c[1])};
  }

  void NelderMead::shrink(double sigma) {
    for (int i = 1; i < 3; ++i) {
      simplex[i].x[0] = simplex[0].x[0] + sigma * (simplex[i].x[0] - simplex[0].x[0]);
      simplex[i].x[1] = simplex[0].x[1] + sigma * (simplex[i].x[1] - simplex[0].x[1]);
      simplex[i].fval = f(simplex[i].x);
    }
  }

}  // namespace pecsim
