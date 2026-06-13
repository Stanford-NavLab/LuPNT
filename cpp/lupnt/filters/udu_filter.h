#pragma once

#include "lupnt/filters/ekf.h"

namespace lupnt {

  /// @brief Extended Kalman filter that stores and updates the covariance as a UDU
  /// (Bierman-Thornton / square-root) factorization `P = U diag(D) U^T` instead of a dense
  /// matrix.
  ///
  /// A drop-in `EKF` replacement (registered with the `Filter` asset factory as `"UDUEKF"`)
  /// for applications that need the improved numerical conditioning of a square-root
  /// filter; reuses `EKF`'s outlier rejection, fault detection, and smoother
  /// infrastructure, but overrides `Predict`/`Update`/`SetCovariance` to maintain `U_`/
  /// `D_diag_` via `UDUDecomposition`/`UDUReconstruct` (see `udu_utils.h`).
  class UDUEKF : public EKF {
  protected:
    MatXd U_;
    VecXd D_diag_;

    /// @brief Refactor the current dense covariance `P_` into `U_`/`D_diag_` via
    /// `UDUDecomposition`.
    ///
    /// Called by `SetCovariance` (on initialization) and at the end of `Predict` (after
    /// `P_` has been propagated and process noise added).
    void SetUd();

    /// @brief Sequentially process each scalar component of the measurement residual
    /// `dz_` using the Carlson-form (rank-1 UD/Bierman-Thornton) measurement update,
    /// refactoring `U_`/`D_diag_` after each component and accumulating the total state
    /// correction `dx_`.
    ///
    /// Called by `Update` once outlier rejection has finalized `dz_`/`H_`/`R_`. Implements
    /// the algorithmic core that distinguishes `UDUEKF::Update` from `EKF::Update`'s
    /// batch Joseph-form update.
    void CarlsonUpdate();

  public:
    UDUEKF() = default;

    /// @brief Construct a UDU EKF from a YAML config node (forwarded to `EKF`).
    explicit UDUEKF(Config& config) : EKF(config) {}
    ~UDUEKF() override = default;

    /// @brief Set the initial state covariance `P_` (via `Filter::SetCovariance`) and
    /// immediately factor it into `U_`/`D_diag_` via `SetUd`.
    void SetCovariance(const MatXd& P) override;

    /// @brief UDU-filter predict step: same propagation as `EKF::Predict` (state via
    /// `f_dyn_`/`F_`, process noise `Q_` via `f_proc_`), but reconstructs `P_` from
    /// `U_`/`D_diag_` before propagation and refactors the result back into `U_`/`D_diag_`
    /// via `SetUd` afterward.
    void Predict(Real t, const State* u = nullptr) override;

    /// @brief UDU-filter update step: same outlier handling as `EKF::Update`, but applies
    /// the correction via the sequential Carlson (UD/Bierman-Thornton) update
    /// (`CarlsonUpdate`) instead of a batch Joseph-form Kalman update.
    void Update(const VecX& z_true) override;

    /// @brief Get the unit-upper-triangular `U` factor of the covariance, `P_ = U diag(D)
    /// U^T`, size [n_x x n_x].
    MatXd GetUFactor() const { return U_; }

    /// @brief Get the diagonal `D` factor of the covariance, `P_ = U diag(D) U^T`,
    /// size [n_x].
    VecXd GetDFactor() const { return D_diag_; }
  };

  /// @brief UDU EKF with one-step stochastic cloning for time-differenced measurements.
  ///
  /// The filter state is `[x_current; x_previous]`. During `Predict`, the first block is
  /// propagated from `t_{k-1}` to `t_k`, while the second block is set to the previous
  /// posterior current state. This lets measurement functions build rows that depend on both
  /// the current and previous receiver state, such as TDCP.
  class UDUStochasticCloningEKF : public UDUEKF {
  public:
    UDUStochasticCloningEKF() = default;
    explicit UDUStochasticCloningEKF(Config& config) : UDUEKF(config) {}
    ~UDUStochasticCloningEKF() override = default;

    void SetBaseStateSize(int base_state_size);
    int GetBaseStateSize() const { return base_state_size_; }

    void SetState(const State& x) override;
    void SetCovariance(const MatXd& P) override;
    void Predict(Real t, const State* u = nullptr) override;

    State GetCurrentState() const { return x_.head(base_state_size_); }
    MatXd GetCurrentCovariance() const {
      return P_.topLeftCorner(base_state_size_, base_state_size_);
    }

  private:
    int base_state_size_ = 0;

    void InferBaseStateSize();
  };

}  // namespace lupnt
