#pragma once

#include <Eigen/Dense>

#include "lupnt/numerics/filters/ekf.h"

namespace lupnt {

  /// @brief Square-Root Information Filter (SRIF): an `EKF` that stores and updates
  /// the covariance as an upper-triangular information square root `R`
  /// (`R^T R = P^{-1}`, so `P = R^{-1} R^{-T}`) instead of a dense matrix.
  ///
  /// A drop-in `EKF` replacement (registered with the `Filter` asset factory as
  /// `"SRIF"`), analogous to `UDUEKF` but factoring the *information* matrix rather
  /// than the covariance. `Predict` and `Update` are performed by orthogonal
  /// (Householder QR) triangularizations of the information array, which are far
  /// better conditioned than the covariance-form Kalman recursions -- the estimate
  /// stays consistent even where the geometry is nearly unobservable (e.g. an
  /// unobserved coast between tracking passes). The dense `P_`/`x_` maintained by the
  /// base class are kept in sync after every step (`P_ = R^{-1} R^{-T}`), so the
  /// inherited outlier rejection, fault detection, logging, and Rauch-Tung-Striebel
  /// (RTS) smoother (`InitializeLogger`/`InitializeSmootherState`/`LogFilterEstimate`/
  /// `UpdateSmoother`) all work unchanged. The filter re-linearizes about its own
  /// estimate each step (extended SRIF), like the EKF.
  class SRIF : public EKF {
  protected:
    /// Upper-triangular information square root, `R_info_^T R_info_ = P_^{-1}`.
    MatXd R_info_;

    /// @brief Factor the current dense covariance `P_` into `R_info_` via a Cholesky
    /// of `P_^{-1}` (so `R_info_` is upper-triangular with `R_info_^T R_info_ = P_^{-1}`).
    ///
    /// Called by `SetCovariance` (on initialization) and is kept consistent with `P_`
    /// by `Predict`/`Update` thereafter.
    void SetInfoSqrt();

    /// @brief Information-form measurement update on the finalized `dz_`/`H_`/`R_`:
    /// whitens the rows by `R_^{-1/2}`, folds them into `R_info_` by a QR
    /// triangularization of `[[R_info_]; [H_w]]`, applies the resulting state
    /// correction `dx_` to `x_`, and reconstructs `P_` from the updated `R_info_`.
    ///
    /// Called by `Update` once outlier rejection has finalized `dz_`/`H_`/`R_`.
    void InformationUpdate();

  public:
    SRIF() = default;

    /// @brief Construct an SRIF from a YAML config node (forwarded to `EKF`).
    explicit SRIF(Config& config) : EKF(config) {}
    ~SRIF() override = default;

    /// @brief Set the initial state covariance `P_` (via `Filter::SetCovariance`) and
    /// immediately factor it into the information square root `R_info_`.
    void SetCovariance(const MatXd& P) override;

    /// @brief SRIF predict step: same propagation as `EKF::Predict` (state via
    /// `f_dyn_`/`F_`, process noise `Q_` via `f_proc_`), but the a-priori information
    /// square root is formed by a Householder QR of the Dyer-McReynolds time-update
    /// array instead of forming `F P F^T + Q` directly. `P_` is reconstructed from
    /// `R_info_` afterward (equal to `F P F^T + Q` for a full-rank `Q`).
    void Predict(Real t, const State* u = nullptr) override;

    /// @brief SRIF update step: same outlier handling as `EKF::Update`, but the
    /// correction is applied via the square-root information measurement update
    /// (`InformationUpdate`) instead of a batch Joseph-form Kalman update.
    void Update(const VecX& z_true) override;

    /// @brief Get the upper-triangular information square root `R` of the covariance,
    /// `P_ = R^{-1} R^{-T}`, size [n_x x n_x].
    MatXd GetInfoSqrt() const { return R_info_; }
  };

  /// @brief Discrete continuous-white-noise-acceleration (CWNA) process-noise
  /// covariance for a Cartesian `[r (3); v (3)]` state over a step `dt`, with
  /// per-axis acceleration power spectral density `accel_psd` [m^2/s^3]:
  ///
  ///   Q = accel_psd * [[dt^3/3 I3, dt^2/2 I3], [dt^2/2 I3, dt I3]].
  ///
  /// A convenient `ProcessNoiseFunction` payload for orbit-determination filters that
  /// want to account for unmodeled dynamics (higher-order gravity, SRP mismodeling,
  /// ...). Returns a 6x6 matrix.
  MatXd CwnaProcessNoise(double dt, double accel_psd);

}  // namespace lupnt
