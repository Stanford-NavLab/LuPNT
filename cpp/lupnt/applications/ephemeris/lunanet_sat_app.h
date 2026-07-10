#pragma once

#include <string>
#include <utility>
#include <vector>

#include "lupnt/applications/application.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  class LunaNetSatApp;

  /// @brief Interface for modular LunaNet satellite application functions.
  ///
  /// A sub-app owns one piece of spacecraft autonomy or navigation logic, such
  /// as ODTS, orbit prediction, ephemeris fitting, or navigation-message
  /// generation. `LunaNetSatApp` calls each registered sub-app from the main
  /// LuPNT simulation schedule.
  class LunaNetSubApp {
  public:
    LunaNetSubApp() = default;
    explicit LunaNetSubApp(std::string name) : name_(std::move(name)) {}
    virtual ~LunaNetSubApp() = default;

    virtual void Setup(LunaNetSatApp& app);
    virtual void Step(Real t) = 0;
    virtual void Finish();

    const std::string& GetName() const { return name_; }
    void SetName(std::string name) { name_ = std::move(name); }

  protected:
    std::string name_ = "lunanet_subapp";
  };

  /// @brief Composite application attached to a LunaNet satellite agent.
  ///
  /// The simulation schedules this application like any other `Application`.
  /// Each scheduled call fans out to registered sub-apps in insertion order.
  /// This keeps the satellite-level event cadence in one place while allowing
  /// ODTS, prediction, ephemeris fitting, and nav-message generation to be
  /// developed as separate modules.
  class LunaNetSatApp : public Application {
  public:
    LunaNetSatApp() = default;
    explicit LunaNetSatApp(Config& config);

    void AddSubApp(Ptr<LunaNetSubApp> app);
    const std::vector<Ptr<LunaNetSubApp>>& GetSubApps() const { return sub_apps_; }

    void Setup() override;
    void Step(Real t) override;
    void Finish();

  private:
    std::vector<Ptr<LunaNetSubApp>> sub_apps_;
  };

}  // namespace lupnt
