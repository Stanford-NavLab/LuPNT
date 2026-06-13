#include "src/udu_utils.h"

namespace filtering_sim {
  using namespace lupnt;

  VecMatPair UDUDecomposition(const MatXd& P_in) {
    const int n = int(P_in.rows());
    VecXd D = VecXd::Zero(n);
    MatXd U = MatXd::Zero(n, n);
    MatXd P = P_in;

    const double eps = 1e-16;

    // std::cout << "[UDU Decomposition]" << std::endl;
    // std::cout << "P_in: " << std::endl << P_in << std::endl;

    // Last column (j = n-1)
    D(n - 1) = P(n - 1, n - 1);
    U(n - 1, n - 1) = 1.0;

    if (std::abs(D(n - 1)) > eps) {
      for (int i = 0; i < n - 1; ++i) {
        U(i, n - 1) = P(i, n - 1) / D(n - 1);
      }
    } else {
      D(n - 1) = 0.0;
      // U(i,n-1) can stay 0 since it multiplies D(n-1)=0 anyway
    }

    // Remaining columns j = n-2 ... 0
    for (int j = n - 2; j >= 0; --j) {
      // D(j) = P(j,j) - sum_{k=j+1}^{n-1} U(j,k)^2 D(k)
      double sum = 0.0;
      for (int k = j + 1; k < n; ++k) {
        sum += U(j, k) * D(k) * U(j, k);
      }
      D(j) = P(j, j) - sum;

      U(j, j) = 1.0;

      // U(i,j) for i = j-1 ... 0
      if (std::abs(D(j)) > eps) {
        for (int i = j - 1; i >= 0; --i) {
          double s = 0.0;
          for (int k = j + 1; k < n; ++k) {
            s += U(i, k) * D(k) * U(j, k);
          }
          U(i, j) = (P(i, j) - s) / D(j);  // <-- FIXED
        }
      } else {
        D(j) = 0.0;
        // U(i,j) can stay 0 because it multiplies D(j)=0 in UDU^T
      }
    }

    // std::cout << "P_out: " << std::endl << U * D.asDiagonal() * U.transpose() << std::endl;

    return {D, U};
  }

  VecMatPair ModifiedGramSchmidt(const MatXd& D_tilde, const MatXd& Y) {
    int n = int(Y.rows());  // size of the state vector
    int n_plus_m = int(Y.cols());
    int m = n_plus_m - n;  // size of the process noise (in our filter n=m)

    std::vector<VecXd> b(n + m);
    VecXd D_m(n);
    MatXd U_m = MatXd::Identity(n, n);

    for (int i = 0; i < n; i++) {
      int l = (n - 1) - i;
      b[l] = Y.row(l);
    }

    for (int k = 0; k < n; k++) {
      int j = (n - 1) - k;
      VecXd fj = D_tilde * b[j];
      D_m[j] = b[j].dot(fj);
      fj = fj / D_m[j];

      for (int i = 0; i < j; i++) {
        U_m(i, j) = b[i].dot(fj);
        b[i] = b[i] - U_m(i, j) * b[j];
      }

      U_m(0, 0) = 1.0;
      VecXd f0 = D_tilde * b[0];
      D_m[0] = b[0].dot(f0);
    }

    return {D_m, U_m};
  }

  // VecMatPair AgeeTurnerRankOneUpdate(const MatXd& U_prev, const VecXd& D_prev, const VecXd& a_in,
  //                                    double c) {
  //   // Agee-Turner rank-1 update for UDU decomposition
  //   // Udiag(D)U^T = U_prev * diag(D_prev) * U_prev^T + c * a * a^T
  //   int n = int(D_prev.size());
  //   VecXd C = VecXd::Zero(n);
  //   C(n - 1) = c;

  //   MatXd U_new = MatXd::Zero(n, n);
  //   VecXd D_new = VecXd::Zero(n);

  //   VecXd a = a_in;

  //   for (int j = n - 1; j >= 1; j--) {
  //     U_new(j, j) = 1.0;
  //     D_new(j) = D_prev(j) + C(j) * a(j) * a(j);
  //     for (int k = 0; k < j; k++) {
  //       a(k) = a(k) - U_prev(k, j) * a(j);
  //       U_new(k, j) = U_prev(k, j) + (C(j) * a(k) * a(j)) / D_new(j);
  //     }
  //     C(j - 1) = C(j) * D_prev(j) / D_new(j);
  //   }

  //   D_new(0) = D_prev(0) + C(0) * a(0) * a(0);

  //   return {D_new, U_new};
  // }

  VecMatPair AgeeTurnerRankOneUpdate(const MatXd& U_prev, const VecXd& D_prev, const VecXd& a_in,
                                     double c) {
    int n = int(D_prev.size());
    MatXd U_new = U_prev;
    VecXd D_new = D_prev;

    // Transform update vector
    VecXd w = U_prev.transpose() * a_in;

    // Use a small epsilon relative to double precision
    // 1e-16 is roughly machine epsilon for doubles.
    // We use a slightly larger buffer to avoid denormalized numbers.
    const double EPS = 1e-15;

    for (int j = 0; j < n; j++) {
      double w_j = w(j);
      double d_old = D_new(j);

      // 1. Compute new diagonal
      double d_new = d_old + c * w_j * w_j;

      // 2. Stability Check
      // If d_new is effectively zero, we cannot divide by it.
      // This implies this state component has zero uncertainty.
      if (std::abs(d_new) < EPS) {
        // Case A: Variance collapsed or was already 0.
        // We clamp it to 0 (or EPS) to avoid NaN.
        D_new(j) = 0.0;

        // If the diagonal is 0, this dimension cannot absorb the update.
        // We stop propagating 'c' (energy) down this path to avoid blowing up U.
        // beta = 0, c_next = 0
        c = 0.0;

        // Note: If d_new is negative, your covariance is no longer
        // Positive Definite. In strict filters, you might want to throw or reset:
        // if (d_new < -EPS) throw std::runtime_error("Filter Divergence: Negative D");
      } else {
        // Case B: Normal update
        D_new(j) = d_new;

        // Calculate beta and next c
        double beta = (c * w_j) / d_new;
        double c_next = (c * d_old) / d_new;

        // Update U columns
        for (int k = j + 1; k < n; k++) {
          U_new(j, k) += beta * w(k);
        }

        // Propagate c
        c = c_next;
      }
    }

    return {D_new, U_new};
  }

  MatXd BackwardSubstitution(const MatXd& U, const MatXd& B) {
    int n = int(U.rows());
    int m = int(B.cols());

    MatXd X = B;
    for (int k = n - 1; k >= 0; k--) {
      for (int j = k + 1; j < n; j++) {
        // X.row(k) -= U(k, j) * X.row(j);
        for (int col = 0; col < m; col++) {
          X(k, col) -= U(k, j) * X(j, col);
        }
      }
    }

    return X;
  }

  MatXd ForwardSubstitution(const MatXd& U, const MatXd& B) {
    int n = int(U.rows());
    int m = int(B.cols());
    MatXd X = B;
    for (int k = 0; k < n; k++) {
      for (int j = 0; j < k; j++) {
        // X.row(k) -= U(j, k) * X.row(j); // Note: U(j,k) is (U^T)(k,j)
        for (int col = 0; col < m; col++) {
          X(k, col) -= U(j, k) * X(j, col);
        }
      }
    }

    return X;
  }

  // Solves (U * D * U^T) * X = B for X
  // U: Upper unit triangular (n x n)
  // D: Diagonal (n)
  // B: RHS matrix (n x m)
  MatXd UDUSolve(const MatXd& U, const VecXd& D, const MatXd& B) {
    int n = int(D.size());
    int m = int(B.cols());

    // System: U * D * U^T * X = B
    // Let Y = D * U^T * X
    // Step 1: Solve U * Y = B (Backward Substitution)
    MatXd Y = BackwardSubstitution(U, B);

    // Step 2: Solve D * Z = Y (Diagonal Scaling)
    // Z = U^T * X
    MatXd Z = Y;
    for (int k = 0; k < n; k++) {
      for (int col = 0; col < m; col++) {
        Z(k, col) /= D(k);
      }
    }

    // Step 3: Solve U^T * X = Z (Forward Substitution)
    // U^T is Lower Unit Triangular
    MatXd X = ForwardSubstitution(U, Z);

    return X;
  }
}  // namespace filtering_sim
