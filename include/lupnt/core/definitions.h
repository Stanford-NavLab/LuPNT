/**
 * @file constants.h
 * @author Stanford NAV LAB
 * @brief List of constants
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#define ASSERT_WITH_MESSAGE(condition, message) \
  if (!(condition)) {                           \
    std::ostringstream oss;                     \
    oss << message;                             \
    throw std::runtime_error(oss.str());        \
  }

#define DEFINE_STATIC_VECTORS_MATRICES(size)       \
  using Vec##size = Matrix<Real, size, 1>;         \
  using Vec##size##d = Matrix<double, size, 1>;    \
  using Vec##size##i = Matrix<int, size, 1>;       \
  using Mat##size = Matrix<Real, size, size>;      \
  using Mat##size##d = Matrix<double, size, size>; \
  using Mat##size##i = Matrix<int, size, size>;    \
  using RowVec##size##d = Matrix<double, 1, size>; \
  using RowVec##size##i = Matrix<int, 1, size>;    \
  using RowVec##size = Matrix<Real, 1, size>;

#define DEFINE_DYNAMIC_VECTORS_MATRICES()                       \
  using VecX = Matrix<Real, Eigen::Dynamic, 1>;                 \
  using VecXd = Matrix<double, Eigen::Dynamic, 1>;              \
  using VecXi = Matrix<int, Eigen::Dynamic, 1>;                 \
  using MatX = Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;    \
  using MatXd = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; \
  using MatXi = Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;    \
  using MatX6 = Matrix<Real, Eigen::Dynamic, 6>;                \
  using MatX6d = Matrix<double, Eigen::Dynamic, 6>;             \
  using MatX3 = Matrix<Real, Eigen::Dynamic, 3>;                \
  using MatX3d = Matrix<double, Eigen::Dynamic, 3>;             \
  using RowVecX = Matrix<Real, 1, Eigen::Dynamic>;              \
  using RowVecXd = Matrix<double, 1, Eigen::Dynamic>;           \
  using RowVecXi = Matrix<int, 1, Eigen::Dynamic>;

#define DEFINE_VECTORS_MATRICES()    \
  DEFINE_STATIC_VECTORS_MATRICES(1)  \
  DEFINE_STATIC_VECTORS_MATRICES(2)  \
  DEFINE_STATIC_VECTORS_MATRICES(3)  \
  DEFINE_STATIC_VECTORS_MATRICES(4)  \
  DEFINE_STATIC_VECTORS_MATRICES(5)  \
  DEFINE_STATIC_VECTORS_MATRICES(6)  \
  DEFINE_STATIC_VECTORS_MATRICES(7)  \
  DEFINE_STATIC_VECTORS_MATRICES(8)  \
  DEFINE_STATIC_VECTORS_MATRICES(9)  \
  DEFINE_STATIC_VECTORS_MATRICES(10) \
  DEFINE_DYNAMIC_VECTORS_MATRICES()

namespace lupnt {

  template <typename T> using Ptr = std::shared_ptr<T>;
  template <typename T, typename... Args> Ptr<T> MakePtr(Args&&... args) {
    return std::make_shared<T>(std::forward<Args>(args)...);
  }

  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixX;
  using Eigen::Vector;
  using Eigen::VectorX;

  using Real = autodiff::real;
  template <int rows, int cols> using Mat = Matrix<Real, rows, cols>;
  template <int rows, int cols> using Matd = Matrix<double, rows, cols>;
  template <int size> using Vec = Matrix<Real, size, 1>;
  template <int size> using Vecd = Matrix<double, size, 1>;
  template <int size> using RowVecd = Matrix<double, 1, size>;
  using RowVecXd = Matrix<double, 1, Eigen::Dynamic>;
  using Quat = Eigen::Quaternion<Real>;
  using Quatd = Eigen::Quaternion<double>;
  using AngleAxis = Eigen::AngleAxis<Real>;
  using AngleAxisd = Eigen::AngleAxis<double>;

  DEFINE_VECTORS_MATRICES()

  static Eigen::IOFormat FMT_CLEAN(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
  static Eigen::IOFormat FMT_HEAVY(Eigen::FullPrecision, 0, ", ", ",\n", "[", "]", "[", "]");
  static Eigen::IOFormat FMT_COMPACT(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ";\n ", "",
                                     "", "[", "]");

}  // namespace lupnt
