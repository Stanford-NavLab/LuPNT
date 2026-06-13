#include "lupnt/applications/lunanet_sat_app.h"

#include <fmt/format.h>

#include "lupnt/core/asset_factory.h"
#include "lupnt/core/logger.h"

namespace lupnt {

  void LunaNetSubApp::Setup(LunaNetSatApp& app) {
    Logger::Debug(fmt::format("Setting up {} under {}", name_, app.GetName()), "LunaNetSubApp");
  }

  void LunaNetSubApp::Finish() {
    Logger::Debug(fmt::format("Finishing {}", name_), "LunaNetSubApp");
  }

  LunaNetSatApp::LunaNetSatApp(Config& config) : Application(config) {
    Logger::Debug(fmt::format("Creating {}", name_), "LunaNetSatApp");
  }

  void LunaNetSatApp::AddSubApp(Ptr<LunaNetSubApp> app) {
    LUPNT_CHECK(app != nullptr, "Cannot add null LunaNet sub-app", "LunaNetSatApp");
    sub_apps_.push_back(std::move(app));
  }

  void LunaNetSatApp::Setup() {
    Logger::Debug(fmt::format("Setting up {}", name_), "LunaNetSatApp");
    for (auto& app : sub_apps_) app->Setup(*this);
    if (GetAgent() != nullptr) Application::Setup();
  }

  void LunaNetSatApp::Step(Real t) {
    Logger::Debug(fmt::format("Stepping {}", name_), "LunaNetSatApp", t);
    for (auto& app : sub_apps_) app->Step(t);
  }

  void LunaNetSatApp::Finish() {
    Logger::Debug(fmt::format("Finishing {}", name_), "LunaNetSatApp");
    for (auto& app : sub_apps_) app->Finish();
  }

  REGISTER_FACTORY_CLASS(Application, LunaNetSatApp)

}  // namespace lupnt
