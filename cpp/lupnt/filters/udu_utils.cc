#include "lupnt/filters/udu_utils.h"

#include <cmath>

#include "lupnt/core/error.h"

namespace lupnt {

  VecMatPair UDUDecomposition(const MatXd& P_in) {
    LUPNT_CHECK(P_in.rows() == P_in.cols(), "UDU input matrix must be square", "UDU");
    const int n = static_cast<int>(P_in.rows());
    VecXd D = VecXd::Zero(n);
    MatXd U = MatXd::Zero(n, n);
    const double eps = 1.0e-16;

    for (int j = n - 1; j >= 0; --j) {
      double sum = 0.0;
      for (int k = j + 1; k < n; ++k) {
        sum += U(j, k) * D(k) * U(j, k);
      }
      D(j) = P_in(j, j) - sum;
      U(j, j) = 1.0;

      if (std::abs(D(j)) <= eps) {
        D(j) = 0.0;
        continue;
      }

      for (int i = j - 1; i >= 0; --i) {
        double offdiag_sum = 0.0;
        for (int k = j + 1; k < n; ++k) {
          offdiag_sum += U(i, k) * D(k) * U(j, k);
        }
        U(i, j) = (P_in(i, j) - offdiag_sum) / D(j);
      }
    }

    return {D, U};
  }

  MatXd UDUReconstruct(const MatXd& U, const VecXd& D) {
    LUPNT_CHECK(U.rows() == U.cols(), "UDU U matrix must be square", "UDU");
    LUPNT_CHECK(U.rows() == D.size(), "UDU factor dimensions do not match", "UDU");
    return U * D.asDiagonal() * U.transpose();
  }

  namespace {
    MatXd BackwardSubstitutionUnitUpper(const MatXd& U, const MatXd& B) {
      const int n = static_cast<int>(U.rows());
      MatXd X = B;
      for (int row = n - 1; row >= 0; --row) {
        for (int k = row + 1; k < n; ++k) {
          X.row(row) -= U(row, k) * X.row(k);
        }
      }
      return X;
    }

    MatXd ForwardSubstitutionUnitLowerFromUpperTranspose(const MatXd& U, const MatXd& B) {
      const int n = static_cast<int>(U.rows());
      MatXd X = B;
      for (int row = 0; row < n; ++row) {
        for (int k = 0; k < row; ++k) {
          X.row(row) -= U(k, row) * X.row(k);
        }
      }
      return X;
    }
  }  // namespace

  MatXd UDUSolve(const MatXd& U, const VecXd& D, const MatXd& B) {
    LUPNT_CHECK(U.rows() == U.cols(), "UDU U matrix must be square", "UDU");
    LUPNT_CHECK(U.rows() == D.size(), "UDU factor dimensions do not match", "UDU");
    LUPNT_CHECK(U.rows() == B.rows(), "UDU solve RHS dimension does not match", "UDU");

    MatXd Y = BackwardSubstitutionUnitUpper(U, B);
    MatXd Z = Y;
    for (int i = 0; i < D.size(); ++i) {
      LUPNT_CHECK(std::abs(D(i)) > 0.0, "Cannot solve singular UDU system", "UDU");
      Z.row(i) /= D(i);
    }
    return ForwardSubstitutionUnitLowerFromUpperTranspose(U, Z);
  }

}  // namespace lupnt
