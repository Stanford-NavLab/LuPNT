// Python driver for LuPNT's square-root UDU stochastic-cloning EKF
// (`UDUStochasticCloningEKF`, cpp/lupnt/numerics/filters/udu_filter.h).
//
// The filter's dynamics / process-noise / measurement hooks are wired to C++ lambdas that
// read buffers stashed from Python just before each step, so the numerically-stable
// covariance algebra (Thornton MWGS time update + Carlson rank-1 measurement update +
// one-step stochastic cloning for TDCP) runs entirely in C++ with no per-epoch Python
// callback -- while the caller keeps ownership of the (already-C++) state propagation and of
// the measurement geometry.  This mirrors how LunarGnssODTSSimulation drives the same filter.

// lupnt
#include <lupnt/lupnt.h>

// pybind11
#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

namespace {

  /// @brief Python-driven wrapper around `UDUStochasticCloningEKF`.
  ///
  /// State ordering is the caller's; the process-noise mapping is derived from whatever
  /// (block-diagonal) base process-noise covariance the caller supplies, so the wrapper makes
  /// no assumption about which entries are orbit / clock / SRP.
  class ScUdFilter {
  public:
    explicit ScUdFilter(int base_state_size) : n_(base_state_size) {
      ekf_.SetBaseStateSize(n_);
      ekf_.SetOutlierThreshold(1.0e12);  // measurement editing is the caller's responsibility

      // Dynamics: return the base state the caller propagated, with its base STM. The filter
      // passes in the current base estimate (== what the caller propagated), which we ignore.
      ekf_.SetDynamicsFunction([this](const State& x, Real, Real, const State*, MatXd* F) -> State {
        if (F) *F = F_base_;
        State out(x);  // copy frame/labels
        out.head(n_) = x_next_base_.cast<Real>();
        return out;
      });

      // Process noise: constant per fixed step; the mapping G is (re)installed each call, as
      // the base filter reads G_ from the object at predict time.
      ekf_.SetProcessNoiseFunction([this](const State&, Real, Real) -> MatXd {
        ekf_.SetProcessNoiseMappingMatrix(G_);
        return MatXd(Qd_.asDiagonal());
      });

      // Measurement: return the modelled measurement, Jacobian, and (diagonal) noise the
      // caller computed at the current augmented state.
      ekf_.SetMeasurementFunction([this](const State&, MatXd* H, MatXd* R) -> VecXd {
        if (H) *H = H_;
        if (R) *R = R_;
        return z_prior_;
      });
    }

    /// Base process-noise covariance (n x n, block diagonal). UDU-decomposed into the mapping
    /// G <- U and diagonal weight Qd <- D, exactly as the reference filter routes correlated
    /// noise through SetProcessNoiseMappingMatrix.
    void SetProcessNoise(const MatXd& Q_base) {
      VecMatPair du = UDUDecomposition(Q_base);
      Qd_ = du.first;
      G_ = du.second;
    }

    /// Initialise from a base estimate (n) and base covariance (n x n); both blocks of the
    /// clone are set equal (fully correlated), matching MakeClonedState/MakeClonedCovariance.
    void SetInitial(const VecXd& x0, const MatXd& P0, double t0) {
      State xs(n_);
      xs.head(n_) = x0.cast<Real>();
      State x_clone(2 * n_);
      x_clone.head(n_) = xs;
      x_clone.tail(n_) = xs;
      MatXd P_clone(2 * n_, 2 * n_);
      P_clone.topLeftCorner(n_, n_) = P0;
      P_clone.topRightCorner(n_, n_) = P0;
      P_clone.bottomLeftCorner(n_, n_) = P0;
      P_clone.bottomRightCorner(n_, n_) = P0;
      ekf_.SetTime(Real(t0));
      ekf_.SetBaseStateSize(n_);
      ekf_.SetState(x_clone);
      ekf_.SetCovariance(P_clone);
    }

    /// Time update: the caller supplies the propagated base state and base STM (both computed
    /// from the filter's current base estimate, read via `current_state`).
    void Predict(double t, const VecXd& x_next_base, const MatXd& F_base) {
      x_next_base_ = x_next_base;
      F_base_ = F_base;
      ekf_.Predict(Real(t));
    }

    /// Measurement update with the modelled measurement, its Jacobian over the augmented
    /// state (m x 2n), and a diagonal measurement-noise covariance.
    void Update(const VecXd& z_true, const VecXd& z_prior, const MatXd& H, const VecXd& R_diag) {
      z_prior_ = z_prior;
      H_ = H;
      R_ = MatXd(R_diag.asDiagonal());
      ekf_.Update(z_true.cast<Real>());
    }

    VecXd CurrentState() const { return ekf_.GetCurrentState().cast<double>(); }
    MatXd CurrentCovariance() const { return ekf_.GetCurrentCovariance(); }

    // ---- Delayed-state fixed-interval smoothing ------------------------------------------
    /// Allocate the per-epoch history buffers the backward pass reads.
    void InitializeLogger(int max_tidx) { ekf_.InitializeLogger(max_tidx); }
    /// Record the prior/posterior augmented state and covariance at `tidx` (call once per
    /// epoch, after predict/update).
    void LogFilterEstimate(int tidx) { ekf_.LogFilterEstimate(tidx); }
    /// Seed the backward pass with the final base-state posterior.
    void InitializeSmootherState() { ekf_.InitializeSmootherState(); }
    /// One backward delayed-state smoother step (call for tidx = N-2 down to 0).
    void UpdateSmoother(int tidx) { ekf_.UpdateSmoother(tidx); }
    VecXd SmoothedState(int tidx) { return ekf_.GetSmoothedState(tidx); }
    MatXd SmoothedCovariance(int tidx) { return ekf_.GetSmoothedCovariance(tidx); }

  private:
    int n_;
    UDUStochasticCloningEKF ekf_;
    VecXd x_next_base_, Qd_, z_prior_;
    MatXd F_base_, G_, H_, R_;
  };

}  // namespace

void InitScUdFilter(py::module& m) {
  py::class_<ScUdFilter>(
      m, "ScUdFilter",
      "Python driver for LuPNT's square-root UDU stochastic-cloning EKF. The caller owns the "
      "state propagation and measurement geometry; the C++ filter performs the Thornton MWGS "
      "time update, Carlson measurement update, and one-step cloning. State ordering is the "
      "caller's; call order per epoch is predict(...) then update(...).")
      .def(py::init<int>(), py::arg("base_state_size"))
      .def("set_process_noise", &ScUdFilter::SetProcessNoise, py::arg("Q_base"),
           "Base (n x n, block-diagonal) process-noise covariance; UDU-decomposed into the "
           "mapping/diagonal the square-root time update needs.")
      .def("set_initial", &ScUdFilter::SetInitial, py::arg("x0"), py::arg("P0"), py::arg("t0"),
           "Initialise the cloned state/covariance from a base estimate x0 (n) and covariance "
           "P0 (n x n).")
      .def("predict", &ScUdFilter::Predict, py::arg("t"), py::arg("x_next_base"), py::arg("F_base"),
           "Time update to t; x_next_base (n) and F_base (n x n) are the base state and STM the "
           "caller propagated from current_state().")
      .def("update", &ScUdFilter::Update, py::arg("z_true"), py::arg("z_prior"), py::arg("H"),
           py::arg("R_diag"),
           "Measurement update: modelled measurement z_prior (m), Jacobian H (m x 2n) over the "
           "augmented [current; previous] state, and diagonal measurement noise R_diag (m).")
      .def("current_state", &ScUdFilter::CurrentState,
           "Current base state estimate (n), the leading block of the augmented state.")
      .def("current_covariance", &ScUdFilter::CurrentCovariance,
           "Current base covariance (n x n), the leading block of the augmented covariance.")
      .def("initialize_logger", &ScUdFilter::InitializeLogger, py::arg("max_tidx"),
           "Allocate smoother history buffers for max_tidx epochs (call before the forward "
           "pass).")
      .def("log_filter_estimate", &ScUdFilter::LogFilterEstimate, py::arg("tidx"),
           "Record the augmented prior/posterior at epoch tidx (call once per epoch after "
           "predict/update).")
      .def("initialize_smoother_state", &ScUdFilter::InitializeSmootherState,
           "Seed the backward pass with the final base-state posterior.")
      .def("update_smoother", &ScUdFilter::UpdateSmoother, py::arg("tidx"),
           "One backward step of the delayed-state (stochastic-cloning) smoother. Unlike RTS, "
           "this uses the measurement-induced cross-covariance the clone carries, so it stays "
           "valid with time-differenced measurements. Call for tidx = N-2 down to 0.")
      .def("smoothed_state", &ScUdFilter::SmoothedState, py::arg("tidx"),
           "Smoothed base state x_{tidx|N} (n).")
      .def("smoothed_covariance", &ScUdFilter::SmoothedCovariance, py::arg("tidx"),
           "Smoothed base covariance P_{tidx|N} (n x n).");
}
