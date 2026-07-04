#pragma once

#include "lupnt/filters/ekf.h"

namespace lupnt {

  /// @brief Schmidt (consider-parameter) Extended Kalman Filter.
  ///
  /// A thin, explicitly-named convenience wrapper around `EKF::SetConsiderStateCount`:
  /// the trailing `n_consider` elements of the state vector are "consider" states,
  /// included in the propagated state/covariance and in the measurement Jacobian
  /// (so their uncertainty -- and its correlation with the estimated states --
  /// is correctly reflected in the covariance), but never corrected by `Update`.
  ///
  /// Typical use: an onboard filter that estimates its own agent's state from a
  /// relative measurement to another agent, while "considering" that other
  /// agent's own (broadcast/prior) state and uncertainty without claiming
  /// estimation authority over it, e.g. `IslOdtsSimulation`
  /// (`lupnt/simulations/IslOdts/isl_odts_simulation.h`).
  class SchmidtEKF : public EKF {
  public:
    SchmidtEKF() = default;

    /// @brief Construct a Schmidt EKF whose trailing `n_consider` state elements
    /// are consider states (see `EKF::SetConsiderStateCount`).
    explicit SchmidtEKF(int n_consider) { SetConsiderStateCount(n_consider); }

    virtual ~SchmidtEKF() = default;
  };

}  // namespace lupnt
