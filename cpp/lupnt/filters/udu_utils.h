#pragma once

#include <utility>

#include "lupnt/core/definitions.h"

namespace lupnt {

  using VecMatPair = std::pair<VecXd, MatXd>;

  /// @brief Factor a symmetric positive-(semi)definite covariance matrix as
  /// `P = U diag(D) U^T`, with `U` unit upper-triangular and `D` diagonal.
  ///
  /// Used by `UDUEKF::SetUd` to (re)factor the filter covariance after initialization and
  /// after each `Predict` step, enabling the square-root-style sequential `CarlsonUpdate`.
  /// Diagonal entries of `D` that are numerically zero (within `1e-16`) are clamped to
  /// exactly zero.
  ///
  /// @param P  Symmetric covariance matrix to factor, size [n x n]
  /// @return   Pair `(D, U)`: diagonal factor `D` (size [n]) and unit-upper-triangular
  ///           factor `U` (size [n x n])
  VecMatPair UDUDecomposition(const MatXd& P);

  /// @brief Reconstruct the covariance matrix `P = U diag(D) U^T` from a UDU factorization.
  ///
  /// Used by `UDUEKF::Predict`/`CarlsonUpdate` to recover the dense covariance from
  /// `U_`/`D_diag_` when needed (e.g. before propagation, or to report `P_`).
  ///
  /// @param U  Unit-upper-triangular factor, size [n x n]
  /// @param D  Diagonal factor, size [n]
  /// @return   Reconstructed covariance matrix `P`, size [n x n]
  MatXd UDUReconstruct(const MatXd& U, const VecXd& D);

  /// @brief Solve the linear system `U diag(D) U^T X = B` for `X`, via unit-triangular
  /// back/forward substitution and a diagonal solve.
  ///
  /// @param U  Unit-upper-triangular factor, size [n x n]
  /// @param D  Diagonal factor, size [n] (must be non-zero in every entry)
  /// @param B  Right-hand side, size [n x m]
  /// @return   Solution `X`, size [n x m]
  MatXd UDUSolve(const MatXd& U, const VecXd& D, const MatXd& B);

}  // namespace lupnt
