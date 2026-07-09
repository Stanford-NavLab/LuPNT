#pragma once

#include <vector>

#include "lupnt/applications/lunanet_sat_app.h"
#include "lupnt/dynamics/joint_orbit_clock_dynamics.h"
#include "lupnt/numerics/filters/schmidt_ekf.h"

namespace lupnt {

  /// @brief Onboard-filter tuning for `IslOdtsApp` (measurement-noise sigmas and the
  /// consider-state layout). These describe how the hub's onboard Schmidt-EKF *weights*
  /// its measurements; the truth/environment side (force models, GPS constellation
  /// geometry, visibility) lives in the driving `IslOdtsSimulation`.
  struct IslOdtsAppParams {
    /// Number of satellites carried in the filter state: `satellites[0]` is the hub
    /// (whose 8-state is estimated), `satellites[1..n_sat-1]` are the linked/consider
    /// satellites (one two-way crosslink and one 8-state consider block each).
    int n_sat = 2;

    // Two-way crosslink range / Doppler (range-rate) measurement noise, applied to
    // every crosslink identically.
    double range_sigma_m = 1.0;
    double range_rate_sigma_mps = 1.0e-3;

    // One-way anchor (surface-station) pseudorange measurement noise, applied to
    // every staged pseudorange (only used at epochs where an anchor is staged).
    double pseudorange_sigma_m = 5.0;

    // Isotropic acceleration process-noise 1-sigma driving each satellite block.
    double process_accel_sigma_mps2 = 1.0e-8;

    // Normalized-residual outlier-rejection threshold (large = effectively disabled).
    double outlier_threshold = 1.0e12;
  };

  /// @brief One epoch of onboard measurements handed to `IslOdtsApp` (via
  /// `StageMeasurements`) before the next `Step()` processes them.
  ///
  /// The two-way crosslink range/range-rate rows are always present (one per link, in
  /// link order). The one-way pseudorange rows are optional aiding: `anchor_pos_mci`
  /// holds the (assumed-known) positions of any known-position anchor transmitters
  /// serving the hub this epoch (e.g. a lunar surface station) in the hub's
  /// Moon-centered inertial frame (`Frame::MOON_CI`), and `anchor_pseudorange_m` their
  /// observed one-way pseudoranges (same length). Leave both empty for a crosslink-only
  /// update.
  struct IslOdtsMeasurementEpoch {
    VecXd crosslink_range_m;         // [n_links]
    VecXd crosslink_range_rate_mps;  // [n_links]

    std::vector<Vec3d> anchor_pos_mci;  // anchor transmitter positions [m], Frame::MOON_CI
    VecXd anchor_pseudorange_m;         // [anchor_pos_mci.size()] observed pseudoranges [m]
  };

  /// @brief Onboard inter-satellite-link (+ optional Earth-GNSS) ODTS flight
  /// application: a reusable `LunaNetSubApp` that owns the hub satellite's Schmidt
  /// Extended Kalman Filter (`SchmidtEKF`, `lupnt/numerics/filters/schmidt_ekf.h`) and runs one
  /// predict/measurement-update cycle per scheduled `Step()`.
  ///
  /// The filter estimates the hub's own 8-state `[r, v, clock_bias, clock_drift]` and
  /// carries each linked satellite's 8-state as a Schmidt "consider" block (propagated
  /// and used in the update, never corrected). Each epoch the driving simulation stages
  /// that epoch's measurement geometry/observations with `StageMeasurements`, then calls
  /// `Step(t)`, which predicts to `t`, rebuilds the (per-epoch-variable) combined
  /// measurement model, and applies the update.
  ///
  /// This is the single onboard *sensor-fusion* point: additional sensors (surface-station
  /// pseudoranges, Doppler beacons, star-tracker/IMU, ...) are fused simply by extending
  /// `IslOdtsMeasurementEpoch` and the combined measurement model, without changing the
  /// hosting `LunaNetSatApp` or the driving simulation's schedule. It mirrors the
  /// `LunarODTSApp` pattern used inside `LunarGnssODTSSimulation`.
  class IslOdtsApp : public LunaNetSubApp {
  public:
    /// State size of one satellite block: `[r(3), v(3), clock_bias, clock_drift]`.
    static constexpr int kSubStateSize = 8;

    IslOdtsApp() : LunaNetSubApp("isl_odts") {}
    explicit IslOdtsApp(IslOdtsAppParams params);

    /// @brief Seed the onboard filter. Must be called before `Setup()`.
    /// @param t0             Initial filter epoch [s, TDB].
    /// @param x0             Initial estimate, `[own(8), consider_1(8), ...]`, size
    ///                       `kSubStateSize * n_sat`.
    /// @param P0             Initial covariance, same layout.
    /// @param filter_dynamics Reduced-order joint orbit-clock dynamics used to
    ///                       propagate the own block and every consider block.
    void Configure(Real t0, const VecXd& x0, const MatXd& P0,
                   Ptr<JointOrbitClockDynamics> filter_dynamics);

    /// @brief Build the Schmidt-EKF (time/state/covariance/dynamics/process-noise) from
    /// the `Configure` inputs. The measurement model is (re)installed per epoch in `Step`.
    void Setup(LunaNetSatApp& app) override;

    /// @brief Predict to `t` [s, TDB] and apply the update using the currently-staged
    /// measurement epoch (which is then consumed).
    void Step(Real t) override;

    void Finish() override;

    /// @brief Stage the next epoch's measurements to be processed on the following `Step`.
    void StageMeasurements(const IslOdtsMeasurementEpoch& meas);

    /// @brief Current onboard estimate `[own(8), consider_*(8)...]`.
    State GetEstimate() const { return filter_->GetState(); }
    /// @brief Current onboard covariance, same layout.
    MatXd GetCovariance() const { return filter_->GetCovariance(); }
    /// @brief Pre-fit measurement residual cached from the most recent `Step`
    /// (`[crosslink range/range-rate..., anchor pseudorange...]`).
    const VecXd& GetPrefitResidual() const { return prefit_resid_; }
    /// @brief The underlying Schmidt-EKF (for advanced inspection).
    const Ptr<SchmidtEKF>& GetFilter() const { return filter_; }

  private:
    IslOdtsAppParams params_;
    Ptr<SchmidtEKF> filter_;
    Ptr<JointOrbitClockDynamics> filter_dynamics_;
    Real t0_ = 0.0;
    VecXd x0_;
    MatXd P0_;
    IslOdtsMeasurementEpoch staged_;
    bool has_staged_ = false;
    VecXd prefit_resid_;
  };

}  // namespace lupnt
