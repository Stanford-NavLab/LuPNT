/**
 * @file nealdermead.h
 * @author Keidai Iiyama
 * @brief Implementation of the Nelder-Mead optimization algorithm in 2D
 * (Used to correct initial angles in ray tracing)
 * @version 0.1
 * @date 2025-02-17
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>

#include "lupnt/environment/plasma/core/definitions.h"  // Vec2d definition

namespace pecsim {

  struct SimplexVertex {
    Vec2d x;
    double fval;
  };

  class NelderMead {
  public:
    NelderMead(std::function<double(Vec2d)> func);
    NelderMead(std::function<double(Vec2d)> func, Vec2d x0, double step = 1.0, bool debug = false);
    Vec2d minimize(int max_iters = 200, double tol = 1e-6, double fval_tol = 0.0,
                   bool debug = false);

    void set_initial_simplex(Vec2d x0, double step = 1.0, bool debug = false);
    void set_initial_simplex(Vec2d x0, double f0, double step = 1.0, bool debug = false);
    int get_iter() const { return iter_; }
    double get_min_val() const { return min_val_; }

  private:
    int iter_ = 0;          // Current iteration number
    double min_val_ = 0.0;  // Minimum function value found
    std::function<double(Vec2d)> f;
    std::array<SimplexVertex, 3> simplex;
    Vec2d mean_point(const std::array<Vec2d, 2>& points);
    Vec2d reflect(const Vec2d& c, const Vec2d& worst, double alpha = 1.0);
    Vec2d expand(const Vec2d& c, const Vec2d& reflected, double gamma = 2.0);
    Vec2d contract(const Vec2d& c, const Vec2d& worst, double rho = 0.5);
    void shrink(double sigma = 0.5);
  };

};  // namespace pecsim
