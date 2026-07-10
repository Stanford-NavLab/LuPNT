#pragma once

#include "lupnt/applications/application.h"
#include "lupnt/core/config.h"
#include "lupnt/simulations/ephemeris/ephemeris_simulation.h"

namespace lupnt {

  /// @brief Parse an `EphemerisSimulationConfig` from an `application:` config block.
  ///
  /// Reads the scalar tuning, the `orbit:` sub-block (Keplerian elements + `coe_frame`),
  /// the frame strings (`propagate_frame`, `output_frame`), and the `fit_window_minutes`
  /// list. Frames are given as their enum names (e.g. `MOON_CI`, `MOON_OP`, `MOON_PA`).
  EphemerisSimulationConfig ConfigToEphemerisConfig(Config& config);

  /// @brief Coordinator application for the ephemeris/almanac datasize-accuracy study
  /// (Example 9), hosted on a thin `EphemerisManager` agent (mirroring ex8's
  /// `SurfaceStationManager` hosting a `GroundOdtsApp`).
  ///
  /// Its single scheduled `Step` runs the whole body of the (former) monolithic
  /// `EphemerisSimulation::Run()`, in the exact same computation order so the numerics are
  /// preserved bit-for-bit (the study is fully deterministic -- no RNG):
  ///   1. propagate a lunar-satellite truth trajectory (Moon gravity + Earth/Sun),
  ///   2. optionally convert it from the inertial propagate frame to the rotating output
  ///      frame (e.g. `MOON_CI` -> `MOON_PA`),
  ///   3. for each `fit_window_minutes` entry, sample `num_windows` windows, fit and
  ///      quantize the `CartesianEphemeris` / `Almanac` models, and size the broadcast bit
  ///      budget.
  ///
  /// The truth dynamics are built app-internally (exactly as the monolith did). Results are
  /// exposed through `GetResults()` (an `EphemerisResults`), identical to the monolith's.
  class EphemerisApp : public Application {
  public:
    EphemerisApp() = default;
    /// @brief Construct from a YAML `application:` block (self-driving, hosted path).
    explicit EphemerisApp(Config& config);
    /// @brief Construct from an `EphemerisSimulationConfig` struct (unit tests / A-B checks).
    explicit EphemerisApp(EphemerisSimulationConfig config);

    /// @brief Schedule the single one-shot `Step` (APPLICATION priority) on the owning
    /// simulation's event queue.
    void Setup() override;

    /// @brief Run the whole datasize/accuracy sweep (lazily initializing on first call, so a
    /// programmatically-set config is honored). Idempotent: subsequent calls are no-ops.
    void Step(Real t) override;

    void Log(Real /*t*/) override {}

    const EphemerisSimulationConfig& GetConfig() const { return cfg_; }
    const EphemerisResults& GetResults() const { return results_; }

  private:
    /// @brief Propagate the truth trajectory and run the fit/quantize sweep, filling
    /// `results_`. Independent of the owning agent (the truth trajectory is built
    /// app-internally). Called once from the first `Step`.
    void Initialize();

    EphemerisSimulationConfig cfg_;
    EphemerisResults results_;
    bool initialized_ = false;
  };

}  // namespace lupnt
