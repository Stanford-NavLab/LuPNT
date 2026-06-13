#pragma once

#include <lupnt/lupnt.h>

namespace filtering_sim {
  using namespace lupnt;

  using VecMatPair = std::pair<VecXd, MatXd>;
  VecMatPair UDUDecomposition(const MatXd& P_in);
  VecMatPair ModifiedGramSchmidt(const MatXd& D_tilde, const MatXd& Y);

  /**
   * Perform Agee-Turner rank-1 update to UDU decomposition
   * @param U_prev Previous U matrix
   * @param D_prev_diag Previous D diagonal
   * @param u_add Column vector to add
   * @param d_add Scalar to add
   * @return Updated D diagonal and U matrix
   */
  VecMatPair AgeeTurnerRankOneUpdate(const MatXd& U_prev, const VecXd& D_prev_diag,
                                     const VecXd& u_add, double d_add);

  /**
   * Solves U * D * U^T * X = B for X
   * @param U: Upper unit triangular (n x n)
   * @param D: Diagonal elements (n)
   * @param B: RHS matrix (n x m)
   * @return Solution matrix X (n x m)
   */
  MatXd UDUSolve(const MatXd& U, const VecXd& D, const MatXd& B);

  /**
   * Solves U * X = B for X
   * @param U: Upper unit triangular (n x n)
   * @param B: RHS matrix (n x m)
   * @return Solution matrix X (n x m)
   */
  MatXd BackwardSubstitution(const MatXd& U, const MatXd& B);

  /**
   * Solves U^T * X = B for X
   * @param U: Upper unit triangular (n x n)
   * @param B: RHS matrix (n x m)
   * @return Solution matrix X (n x m)
   */
  MatXd ForwardSubstitution(const MatXd& U, const MatXd& B);

}  // namespace filtering_sim
