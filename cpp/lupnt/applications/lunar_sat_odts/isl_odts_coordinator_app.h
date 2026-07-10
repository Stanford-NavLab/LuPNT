#pragma once

#include <limits>
#include <string>
#include <tuple>
#include <vector>

#include "lupnt/applications/application.h"
#include "lupnt/core/config.h"
#include "lupnt/dynamics/joint_orbit_clock_dynamics.h"
#include "lupnt/numerics/filters/ekf.h"
#include "lupnt/simulations/isl_odts/isl_odts_simulation.h"

namespace lupnt {

  class IslOdtsApp;
  class LunaNetSatApp;

  /// @brief Parse an `IslOdtsConfig` from an `application:` config block (the inverse of
  /// `IslOdtsConfigToConfig`). Reads the scalar tuning, the `satellites:` list
  /// (name + Cartesian r0/v0 [m, m/s] in `Frame::MOON_CI` + clock bias/drift), and the
  /// optional `surface_stations:` list.
  IslOdtsConfig ConfigToIslOdtsConfig(Config& config);

  /// @brief Serialize an `IslOdtsConfig` into an `application:`-style `Config` block (the
  /// inverse of `ConfigToIslOdtsConfig`), used by `IslOdtsSimulation` to drive the
  /// agent-based path from the struct API.
  Config IslOdtsConfigToConfig(const IslOdtsConfig& cfg);

  /// @brief Coordinator application for the distributed inter-satellite-link (ISL) ODTS
  /// scenario, hosted on a thin `IslOdtsManager` agent (mirroring ex7's
  /// `GroundStationManagerApp` on a `GroundStationManager`).
  ///
  /// This is the single coordination point of the constellation. Its scheduled `Step(t)`
  /// runs the whole per-epoch body of the (former) monolithic `IslOdtsSimulation::Run()`,
  /// in the exact same computation and RNG-draw order so the numerics are preserved
  /// bit-for-bit:
  ///   1. propagate every satellite's truth `JointOrbitClockDynamics` (orbit + clock,
  ///      per-satellite clock-noise seeding),
  ///   2. synthesize the pairwise two-way crosslink range / range-rate / time-transfer /
  ///      frequency-transfer measurements,
  ///   3. synthesize the surface-station one-way pseudorange / Doppler beacon measurements
  ///      (elevation-gated),
  ///   4a. stage each onboard `IslOdtsApp` filter's rows and `Step` it,
  ///   4b. predict/update the centralized ground EKF from all station pseudoranges,
  ///   5. periodically exchange consider-state (Covariance-Intersection or legacy overwrite).
  ///
  /// The truth and reduced-order filter dynamics are built app-internally (the onboard
  /// filters use a reduced gravity field), exactly as the monolith did. Result series are
  /// exposed through `GetResults()` (an `IslOdtsResults`), identical to the monolith's.
  class IslOdtsCoordinatorApp : public Application {
  public:
    IslOdtsCoordinatorApp() = default;
    /// @brief Construct from a YAML `application:` block (self-driving, hosted path).
    explicit IslOdtsCoordinatorApp(Config& config);
    /// @brief Construct from an `IslOdtsConfig` struct (used by the `IslOdtsSimulation`
    /// adapter and unit tests).
    explicit IslOdtsCoordinatorApp(IslOdtsConfig config);

    /// @brief Build the truth/filter/station/centralized machinery, record epoch 0, and
    /// schedule the periodic per-epoch `Step`s (APPLICATION priority) on the owning
    /// simulation's event queue.
    void Setup() override;

    /// @brief Run one per-epoch coordination body at time `t` [s] (epoch index
    /// `k = round(t/dt_s)`).
    void Step(Real t) override;

    void Log(Real /*t*/) override {}

    const IslOdtsConfig& GetConfig() const { return cfg_; }
    const IslOdtsResults& GetResults() const { return results_; }

  private:
    static constexpr int kSubStateSize = 8;

    /// @brief Build all truth/filter/station/centralized state and record epoch 0. Called
    /// once (lazily on the first `Step`, so any programmatic config change between
    /// construction and `run()` is honored).
    void Initialize();
    /// @brief Run the per-epoch body for epoch index `k` (1..N-1).
    void RunEpoch(int k);
    /// @brief Record the onboard filters' estimates/covariances into `results_` at row `k`.
    void RecordFilters(int k);
    /// @brief Record the centralized filter estimate/covariance into `results_` at row `k`.
    void RecordCentral(int k);

    IslOdtsConfig cfg_;
    IslOdtsResults results_;
    bool initialized_ = false;

    // Scenario dimensions / epoch.
    int n_sat_ = 0;
    int n_links_ = 0;
    int n_state_ = 0;
    int n_stn_ = 0;
    int N_ = 0;
    Real t0_tdb_ = 0.0;
    double next_exchange_s_ = std::numeric_limits<double>::infinity();

    // Truth (one JointOrbitClockDynamics per satellite; per-satellite clock seeds).
    std::vector<Ptr<JointOrbitClockDynamics>> dyn_truth_;
    std::vector<State> x_truth_;
    std::vector<VecXd> x_truth0_;

    // Per-filter block ordering: order_[j][b] = global sat index of block b in filter j;
    // pos_[j][g] = block index of global sat g in filter j.
    std::vector<std::vector<int>> order_;
    std::vector<std::vector<int>> pos_;

    // Onboard distributed filters (each an IslOdtsApp hosted in a LunaNetSatApp).
    std::vector<Ptr<IslOdtsApp>> apps_;
    std::vector<Ptr<LunaNetSatApp>> sat_apps_;

    // Surface-station beacons (fixed in MOON_PA).
    std::vector<IslSurfaceStationConfig> stations_;
    std::vector<Vec3> station_bf_;

    // Centralized ground filter.
    Ptr<EKF> central_;
  };

}  // namespace lupnt
