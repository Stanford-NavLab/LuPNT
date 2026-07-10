#pragma once

#include <utility>
#include <vector>

#include "lupnt/core/definitions.h"
#include "lupnt/measurements/measurement.h"
#include "lupnt/measurements/measurement_utils.h"

namespace lupnt {

  /// @brief Two-way inter-satellite crosslink range/range-rate (+ optional one-way
  /// pseudorange to a known-position transmitter) measurement model for a hub satellite.
  ///
  /// The filter state stacks one `[r(3), v(3), clock_bias, clock_drift]` block per
  /// satellite (`Config::sub_state_size`): block 0 is the hub, blocks `1..n_links` are the
  /// linked satellites. The model produces, per epoch:
  ///   - `2 * n_links` clock-independent two-way crosslink rows `[range, range-rate]`
  ///     (`RangeAndRangeRate` between the hub and each linked block), followed by
  ///   - one one-way pseudorange row per known-position anchor transmitter (e.g. a lunar
  ///     surface station with a known clock), `|r_hub - r_anchor| + C * b_hub`.
  ///
  /// The Jacobian is taken by autodiff *directly on the raw `VecX` state* (not on a
  /// reconstructed `State`), matching the deliberate workaround that keeps the derivative
  /// seeds intact for the `SchmidtEKF` update. This model was previously implemented inline
  /// in `IslOdtsApp`; it now lives here so the application only stages data and wires the
  /// model.
  class IslCrosslinkMeasurement : public MeasurementClone<IslCrosslinkMeasurement> {
  public:
    /// @brief Configuration + per-epoch geometry for the combined crosslink/GPS model.
    struct Config {
      int sub_state_size = 8;   ///< per-satellite state block size `[r,v,cb,cd]`
      int n_links = 0;          ///< number of two-way crosslinks (hub <-> linked sat i)
      int idx_position = 0;     ///< hub position offset within its state block
      int idx_velocity = 3;     ///< hub velocity offset
      int idx_clock_bias = 6;   ///< hub clock-bias offset (for the pseudorange rows)
      int idx_clock_drift = 7;  ///< hub clock-drift offset (for the anchor Doppler rows)

      /// Known-position anchor transmitter positions [m], hub Moon-centered inertial
      /// frame (`Frame::MOON_CI`), e.g. a lunar surface station; empty for a
      /// crosslink-only update.
      std::vector<Vec3d> anchor_pos_mci;
      /// Anchor inertial velocities [m/s], `Frame::MOON_CI` (parallel to `anchor_pos_mci`),
      /// required only when `include_anchor_doppler` is set.
      std::vector<Vec3d> anchor_vel_mci;

      double sigma_range_m = 1.0;            ///< two-way crosslink range noise 1-sigma [m]
      double sigma_range_rate_mps = 1.0e-3;  ///< crosslink range-rate noise 1-sigma [m/s]
      double sigma_pseudorange_m = 200.0;    ///< one-way anchor pseudorange noise 1-sigma [m]

      /// When true, each anchor also yields a one-way **Doppler** (pseudorange-rate) row,
      /// `u . (v_hub - v_anchor) + C * clock_drift_hub` (`u` = line of sight), appended
      /// after all anchor pseudorange rows. Adds velocity + clock-drift observability.
      bool include_anchor_doppler = false;
      double sigma_anchor_doppler_mps = 1.0e-3;  ///< one-way anchor Doppler noise 1-sigma [m/s]

      /// When true, each two-way crosslink also yields a **two-way time-transfer** row:
      /// the clock-bias difference `C * (b_hub - b_link_i)` between the two endpoints.
      /// Two-way ranging cancels the clocks (its sum-of-one-ways gives range), whereas the
      /// *difference* of the same two one-way exchanges recovers the relative clock offset.
      /// These rows are appended after the anchor pseudorange rows.
      bool include_time_transfer = false;
      double sigma_time_transfer_m = 1.0;  ///< two-way time-transfer noise 1-sigma [m]

      /// When true, each crosslink also yields a two-way **frequency-transfer** row: the
      /// range-rate-equivalent clock-*drift* difference `C * (d_hub - d_link_i)`, appended
      /// after the time-transfer rows. This is the rate companion to the time transfer and
      /// makes the endpoints' relative clock drift directly observable.
      bool include_frequency_transfer = false;
      double sigma_frequency_transfer_mps = 1.0e-3;  ///< frequency-transfer noise 1-sigma [m/s]
    };

    IslCrosslinkMeasurement() = default;
    explicit IslCrosslinkMeasurement(Config config) : config_(std::move(config)) {}

    const Config& GetConfig() const { return config_; }
    void SetConfig(const Config& config) { config_ = config; }

    MeasData Compute(const State& x, MatXd* H = nullptr) const override;

  private:
    Config config_;

    /// @brief Number of measurement rows `2 * n_links + anchor_pos_mci.size()`.
    int NumRows() const;
    /// @brief Raw (autodiff-able) measurement model `z = h(x)`.
    VecX Model(const VecX& x) const;
    /// @brief Diagonal measurement noise covariance `R`.
    MatXd Covariance() const;
  };

}  // namespace lupnt
