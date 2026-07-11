#pragma once

#include <utility>

#include "lupnt/core/definitions.h"
#include "lupnt/measurements/measurement.h"

namespace lupnt {

  /// @brief Satellite-to-satellite angles-only (bearing) measurement: the unit
  /// line-of-sight direction from an OBSERVER satellite to a TARGET satellite whose
  /// position is treated as KNOWN (its ephemeris is supplied per epoch).
  ///
  /// The filter state is the observer's own orbit 6-state `[r(3), v(3)]` (no clock --
  /// angles are clock-independent). Modelled on `IslCrosslinkMeasurement` for the
  /// `Model(const VecX&)` + autodiff-`Compute` + `Covariance()`/`NumRows()` pattern, and on
  /// `LanderCraterMeasurement` for the unit-line-of-sight residual form. Because the target
  /// position is known, this is a *landmark-style* bearing update: it directly constrains
  /// the observer position (transverse to the line of sight) and, through the orbit dynamics
  /// over the arc, the full 6-state -- unlike the ill-posed problem of estimating the
  /// target's range from angles alone (which is deliberately NOT attempted here).
  ///
  /// `Model(x)` returns the 3-component unit vector `u = (target - r_obs) / |target - r_obs|`
  /// in `Frame::MOON_CI`; the angular noise `sigma_rad` is applied per line-of-sight
  /// component (small-angle unit-vector approximation), so `Covariance() = sigma_rad^2 I_3`.
  class SatBearingMeasurement : public MeasurementClone<SatBearingMeasurement> {
  public:
    /// @brief Configuration + per-epoch geometry for the bearing model.
    struct Config {
      int idx_position = 0;                  ///< observer position offset in the filter state
      Vec3d target_pos_mci = Vec3d::Zero();  ///< known target position [m], Frame::MOON_CI
      double sigma_rad = 1.0e-4;             ///< bearing noise 1-sigma [rad], per LOS component
    };

    SatBearingMeasurement() = default;
    explicit SatBearingMeasurement(Config config) : config_(std::move(config)) {}

    const Config& GetConfig() const { return config_; }
    void SetConfig(const Config& config) { config_ = config; }

    MeasData Compute(const State& x, MatXd* H = nullptr) const override;

  private:
    Config config_;

    /// @brief Number of measurement rows (3: the unit line-of-sight components).
    int NumRows() const { return 3; }
    /// @brief Raw (autodiff-able) measurement model `u = h(x)` (unit line of sight).
    VecX Model(const VecX& x) const;
    /// @brief Diagonal measurement noise covariance `R = sigma_rad^2 I_3`.
    MatXd Covariance() const;
  };

}  // namespace lupnt
