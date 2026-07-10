#pragma once

#include <vector>

#include "lupnt/applications/application.h"
#include "lupnt/core/config.h"
#include "lupnt/simulations/lunar_gnss_odts/lunar_gnss_odts_simulation.h"

namespace lupnt {

  /// @brief Parse a `LunarGnssODTSConfig` from an `application:` config block.
  ///
  /// The block carries the same nested sections as the standalone scenario file
  /// (`simulation:`, `pipeline:`, `truth:`, `constellation:`, `plasma:`, `dynamics:`,
  /// `measurements:`, `filter:`, `receiver_app:`, `design:`); this delegates to
  /// `ParseLunarGnssODTSConfig`, resolving relative data/output paths against the current
  /// working directory (the run script injects absolute paths).
  LunarGnssODTSConfig ConfigToLunarGnssODTSConfig(Config& config);

  /// @brief Coordinator application for the lunar-orbiting GNSS ODTS scenario (Example 6),
  /// hosted on a physical `Spacecraft` receiver agent, which owns the truth orbit+clock and whose
  /// self-propagated truth the engine reads.
  ///
  /// This wraps the (numerics-preserving) free-function ODTS engine. Its scheduled `Step`
  /// runs the entire GNSS Monte-Carlo body once, lazily, in the exact same computation and
  /// RNG-draw order as the former `LunarGnssODTSSimulation::Run()` so the numerics are
  /// preserved bit-for-bit:
  ///   1. build the receiver truth trajectory + cislunar GNSS sidelobe link geometry / CN0
  ///      (from the precompute cache when available, else in-memory), applying the optional
  ///      plasma delay table,
  ///   2. run `monte_carlo_runs` seeds of the UDU EKF (or UDU stochastic-cloning EKF when
  ///      TDCP is enabled), writing `trajectory_mc<N>.csv` + `summary.csv` under
  ///      `output_dir`.
  ///
  /// The heavy engine (`RunLunarGnssODTSMonteCarlo` and its helpers) is unchanged; the app
  /// only resets the global epoch to 0 (so the engine's absolute-TDB propagation is not
  /// double-counted through `GetLupntEpoch()`) and drives it. Per-seed `LunarGnssODTSSummary`
  /// series are exposed via `GetSummaries()`; the full time series live in the CSV outputs.
  class LunarGnssOdtsApp : public Application {
  public:
    LunarGnssOdtsApp() = default;
    /// @brief Construct from a YAML `application:` block (self-driving, hosted path).
    explicit LunarGnssOdtsApp(Config& config);
    /// @brief Construct from a `LunarGnssODTSConfig` struct (struct API / tests).
    explicit LunarGnssOdtsApp(LunarGnssODTSConfig config);

    /// @brief Schedule the single per-run `Step` (APPLICATION priority) on the owning
    /// simulation's event queue. The heavy engine is deferred to that `Step`.
    void Setup() override;

    /// @brief Run the whole Monte-Carlo ODTS body once (lazy; subsequent calls are no-ops).
    void Step(Real t) override;

    void Log(Real /*t*/) override {}

    /// @brief Build the receiver truth trajectory + link geometry cache without running the
    /// EKF, writing `links_file`. Optional; the EKF `Step` computes links in-memory when no
    /// valid cache exists.
    void Precompute();

    const LunarGnssODTSConfig& GetConfig() const { return cfg_; }
    const std::vector<LunarGnssODTSSummary>& GetSummaries() const { return summaries_; }

  private:
    /// @brief Reset the global epoch and drive `RunLunarGnssODTSMonteCarlo(cfg_)`.
    void RunAll();

    LunarGnssODTSConfig cfg_;
    std::vector<LunarGnssODTSSummary> summaries_;
    bool ran_ = false;
  };

}  // namespace lupnt
