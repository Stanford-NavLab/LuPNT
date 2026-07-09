#include "lupnt/numerics/filters/udu_utils.h"

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

  VecMatPair ModifiedGramSchmidt(const MatXd& D_tilde, const MatXd& Y) {
    const int n = static_cast<int>(Y.rows());  // state dimension
    LUPNT_CHECK(D_tilde.rows() == Y.cols() && D_tilde.cols() == Y.cols(),
                "MWGS weight matrix must be square with size Y.cols()", "UDU");

    std::vector<VecXd> b(n);
    for (int i = 0; i < n; ++i) b[i] = Y.row(i).transpose();

    VecXd D_m = VecXd::Zero(n);
    MatXd U_m = MatXd::Identity(n, n);

    // Orthogonalize from the last row toward the first (Thornton ordering). Unlike
    // reconstructing P and re-factoring, this never forms the full covariance, so diagonal
    // entries many orders of magnitude below the largest are preserved rather than rounded
    // away. Only exactly-zero/negative directions are clamped (a tiny positive variance --
    // e.g. a clock-drift term -- must survive).
    for (int k = 0; k < n; ++k) {
      const int j = (n - 1) - k;
      VecXd fj = D_tilde * b[j];
      D_m(j) = b[j].dot(fj);
      if (D_m(j) <= 0.0) {
        D_m(j) = 0.0;
        continue;
      }
      fj /= D_m(j);
      for (int i = 0; i < j; ++i) {
        U_m(i, j) = b[i].dot(fj);
        b[i] -= U_m(i, j) * b[j];
      }
    }
    return {D_m, U_m};
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
