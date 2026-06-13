#pragma once

#include "lupnt/applications/application.h"
#include "lupnt/core/config.h"
#include "lupnt/filters/filter.h"

namespace lupnt {

  /// @brief Application that runs a lunar/Earth surface ground station's tracking and
  /// navigation-filter logic.
  ///
  /// Currently a minimal skeleton (Setup()/Step() only log debug messages); intended to be
  /// extended to drive `filter_` over the station's tracking measurements, analogous to how
  /// RoverApp drives its filter.
  class SurfaceStationApp : public Application {
  protected:
    Ptr<Filter> filter_;

  public:
    SurfaceStationApp() = default;

    /// @brief Construct a SurfaceStationApp from a YAML Config node.
    ///
    /// Called via AssetFactory<Application, Config&>::Create from the Agent(Config&)
    /// constructor when the agent's `application:` block has `class: SurfaceStationApp`.
    ///
    /// @param config YAML configuration node (forwarded to Application(Config&), e.g.
    ///                `name`/`frequency`).
    SurfaceStationApp(Config& config);

    /// @brief SurfaceStationApp-specific override of Application::Step(): currently a
    /// placeholder that only logs a debug message at time `t`.
    ///
    /// Invoked once per scheduling period (set up in Setup()) by the Simulation event loop.
    ///
    /// @param t Current simulation time [s, since simulation epoch]
    void Step(Real t) override;

    /// @brief SurfaceStationApp-specific override of Application::Setup(): currently a
    /// placeholder that only logs a debug message; does not yet call the base Setup() or
    /// initialize `filter_`.
    void Setup() override;
  };

}  // namespace lupnt
