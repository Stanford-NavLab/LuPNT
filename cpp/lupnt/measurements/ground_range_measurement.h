#pragma once

#include <utility>

#include "lupnt/core/definitions.h"
#include "lupnt/measurements/measurement.h"
#include "lupnt/measurements/measurement_utils.h"

namespace lupnt {

  /// @brief Instantaneous range and/or range-rate of a target state relative to a fixed
  /// reference state (e.g. a ground station), with the closed-form design matrix.
  ///
  /// Given a target `[r(3), v(3)]` (the `State` passed to `Compute`) and a reference
  /// `[r(3), v(3)]` (`Config::reference_state`, same frame), it returns the selected
  /// observables (`range = |dr|`, `range-rate = dr.dv/|dr|`, with `dr = r - r_ref`) and, when
  /// a Jacobian is requested, the instantaneous partials with respect to the target Cart6:
  ///   - range row     : `[u^T, 0]` with `u = dr/|dr|`
  ///   - range-rate row : `[((dv - rdot*u)/rho)^T, u^T]`.
  ///
  /// Batch orbit-determination chains this instantaneous design matrix with the state
  /// transition matrix `Phi(t_k, t0)` outside the model; this class owns only the geometry.
  class GroundStationRangeMeasurement : public MeasurementClone<GroundStationRangeMeasurement> {
  public:
    struct Config {
      bool use_range = true;       ///< include the range row
      bool use_range_rate = true;  ///< include the range-rate row
      /// Reference (station) Cart6 `[r, v]` in the same frame as the target state.
      Vec6d reference_state = Vec6d::Zero();
      double range_sigma_m = 1.0;            ///< range noise 1-sigma [m]
      double range_rate_sigma_mps = 1.0e-3;  ///< range-rate noise 1-sigma [m/s]
    };

    GroundStationRangeMeasurement() = default;
    explicit GroundStationRangeMeasurement(Config config) : config_(std::move(config)) {}

    const Config& GetConfig() const { return config_; }
    void SetConfig(const Config& config) { config_ = config; }

    MeasData Compute(const State& x, MatXd* H = nullptr) const override;

  private:
    Config config_;
    int NumRows() const { return (config_.use_range ? 1 : 0) + (config_.use_range_rate ? 1 : 0); }
  };

}  // namespace lupnt
