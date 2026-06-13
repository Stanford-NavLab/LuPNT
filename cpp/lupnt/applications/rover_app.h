#pragma once

#include "lupnt/applications/application.h"
#include "lupnt/core/config.h"
#include "lupnt/dynamics/dynamics.h"
#include "lupnt/filters/filter.h"

namespace lupnt {

  /// @brief Application that drives a 2D ground rover's navigation filter (e.g. an EKF over
  /// SurfaceDynamics2D) from a configured trajectory command.
  class RoverApp : public Application {
  protected:
    Ptr<Filter> filter_;
    Ptr<Dynamics> dynamics_;

  public:
    RoverApp() = default;

    /// @brief Construct a RoverApp from a YAML Config node.
    ///
    /// Called via AssetFactory<Application, Config&>::Create from the Agent(Config&)
    /// constructor when the agent's `application:` block has `class: RoverApp`. Reads the
    /// required `dynamics` and `filter` sub-configs, prefixes their names with this
    /// application's name, and instantiates the corresponding Dynamics and Filter objects
    /// via their AssetFactory registries. Throws (LUPNT_CHECK) if either sub-config is
    /// missing.
    ///
    /// @param config YAML configuration node with at least `dynamics:` and `filter:` blocks.
    RoverApp(Config& config);

    /// @brief RoverApp-specific override of Application::Step(): computes the rover's
    /// commanded linear/angular velocity from the configured circular trajectory and logs it.
    ///
    /// Invoked once per scheduling period (set up in Setup()) by the Simulation event loop.
    ///
    /// @param t Current simulation time [s, since simulation epoch]
    void Step(Real t) override;

    /// @brief RoverApp-specific override of Application::Setup(): initializes the
    /// navigation filter's initial state, covariance, dynamics-propagation function, and
    /// process-noise function from this rover's configured initial-state uncertainty.
    ///
    /// Builds the initial covariance P0 from `initial_state.sigma_r` [m] (position std.
    /// dev., applied to both x and y) and `initial_state.sigma_theta` [deg] (heading std.
    /// dev., converted to rad), seeds the filter with the agent's current time and state
    /// (Cart6), wires `filter_`'s dynamics-propagation callback to `dynamics_->Propagate`,
    /// and uses P0 as a (constant) process-noise matrix. Called once before the simulation
    /// loop starts.
    void Setup() override;

    /// @brief RoverApp-specific override of Application::Log(): forwards to the navigation
    /// filter's own Log(t) (state/covariance), in addition to the base debug message.
    ///
    /// @param t Current simulation time [s, since simulation epoch]
    void Log(Real t) override;
  };

}  // namespace lupnt
