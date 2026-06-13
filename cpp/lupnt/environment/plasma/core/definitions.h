/**
 * @file eigen_utils.h
 * @author Keidai Iiyama
 * @brief This file contains utility functions for Eigen operations.
 * @version 0.1
 * @date 2025-02-17
 */

#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>

#define DEFINE_VECTORS_MATRICES_FIXED_GCPM(size)              \
  using Vec##size##d = Matrix<double, size, 1>;               \
  using Vec##size##i = Matrix<int, size, 1>;                  \
  using Mat##size##d = Matrix<double, size, size>;            \
  using Mat##size##i = Matrix<int, size, size>;               \
  using RowVec##size##d = Matrix<double, 1, size>;            \
  using RowVec##size##i = Matrix<int, 1, size>;               \
  using MatX##size##d = Matrix<double, Eigen::Dynamic, size>; \
  using Mat##size##Xd = Matrix<double, size, Eigen::Dynamic>;

#define DEFINE_VECTORS_MATRICES_DYNAMIC_GCPM()                  \
  using VecXd = Matrix<double, Eigen::Dynamic, 1>;              \
  using VecXi = Matrix<int, Eigen::Dynamic, 1>;                 \
  using MatXd = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; \
  using MatXi = Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;    \
  using RowVecXd = Matrix<double, 1, Eigen::Dynamic>;           \
  using RowVecXi = Matrix<int, 1, Eigen::Dynamic>;

#define DEFINE_VECTORS_MATRICES_GCPM()   \
  DEFINE_VECTORS_MATRICES_FIXED_GCPM(1)  \
  DEFINE_VECTORS_MATRICES_FIXED_GCPM(2)  \
  DEFINE_VECTORS_MATRICES_FIXED_GCPM(3)  \
  DEFINE_VECTORS_MATRICES_FIXED_GCPM(4)  \
  DEFINE_VECTORS_MATRICES_FIXED_GCPM(5)  \
  DEFINE_VECTORS_MATRICES_FIXED_GCPM(6)  \
  DEFINE_VECTORS_MATRICES_FIXED_GCPM(7)  \
  DEFINE_VECTORS_MATRICES_FIXED_GCPM(8)  \
  DEFINE_VECTORS_MATRICES_FIXED_GCPM(9)  \
  DEFINE_VECTORS_MATRICES_FIXED_GCPM(10) \
  DEFINE_VECTORS_MATRICES_DYNAMIC_GCPM()

namespace pecsim {

  // Eigen
  using Eigen::Block;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixBase;
  using Eigen::MatrixX;
  using Eigen::Vector;
  using Eigen::VectorX;

  template <int rows, int cols> using Matd = Matrix<double, rows, cols>;
  template <int size> using Vecd = Matrix<double, size, 1>;
  template <int size> using RowVecd = Matrix<double, 1, size>;

  // Define fixed-size vectors and matrices for GCPM
  DEFINE_VECTORS_MATRICES_GCPM()

}  // namespace pecsim
