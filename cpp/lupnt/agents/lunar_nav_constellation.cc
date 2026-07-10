#include "lupnt/agents/lunar_nav_constellation.h"

#include "lupnt/core/asset_factory.h"

namespace lupnt {

  void LunarNavConstellation::AddSatellite(const Config& shared, const YAML::Node& initial_state,
                                           const std::string& nm) {
    YAML::Node child = YAML::Clone(shared);  // per-satellite copy of the shared dynamics/clock
    child["name"] = name_ + "/" + nm;
    child["initial_state"] = YAML::Clone(initial_state);
    Config child_cfg(child);
    sats_.push_back(MakePtr<Spacecraft>(child_cfg));
    sat_names_.push_back(nm);
  }

  LunarNavConstellation::LunarNavConstellation(Config& config) : Agent(config) {
    // Shared template: the force model (+ optional clock, seed) every satellite inherits.
    LUPNT_CHECK(config["dynamics"], "LunarNavConstellation needs a shared `dynamics:` block",
                "LunarNavConstellation");
    YAML::Node shared;
    shared["dynamics"] = config["dynamics"];
    if (config["clock"]) shared["clock"] = config["clock"];
    if (config["seed"]) shared["seed"] = config["seed"];

    if (config["satellites"]) {
      // Explicit per-satellite states (each inherits the shared dynamics/clock).
      for (const auto& item : config["satellites"]) {
        YAML::Node s = item;
        std::string nm = s["name"].as<std::string>("SV");
        YAML::Node init;
        if (s["initial_state"]) {
          init = s["initial_state"];
        } else if (s["r0_m"]) {  // {name, r0_m, v0_mps, [clock_bias_s, clock_drift_sps]} shorthand
          init["r0_m"] = s["r0_m"];
          init["v0_mps"] = s["v0_mps"];
          if (s["clock_bias_s"]) init["clock_bias_s"] = s["clock_bias_s"];
          if (s["clock_drift_sps"]) init["clock_drift_sps"] = s["clock_drift_sps"];
        } else {
          LUPNT_CHECK(false, "constellation satellite needs `initial_state` or `r0_m`/`v0_mps`",
                      "LunarNavConstellation");
        }
        AddSatellite(Config(shared), init, nm);
      }
    } else if (config["walker"]) {
      // Symmetric Walker constellation of a common frozen orbit.
      Config w(config["walker"]);
      int n_planes = w["n_planes"].as<int>();
      int per_plane = w["sats_per_plane"].as<int>();
      double a = w["a"].as<double>(), e = w["e"].as<double>(), inc = w["i"].as<double>();
      double omega = w["omega"].as<double>();
      double raan0 = w["raan0_deg"].as<double>(0.0);
      double m0 = w["m0_deg"].as<double>(0.0);
      double phase = w["phase_deg"].as<double>(0.0);  // inter-plane mean-anomaly phasing
      std::string frame = w["frame"].as<std::string>("MOON_OP");
      std::string prefix = w["name_prefix"].as<std::string>("SV");
      int idx = 1;
      for (int p = 0; p < n_planes; ++p) {
        for (int s = 0; s < per_plane; ++s) {
          YAML::Node init;
          init["class"] = "ClassicalOE";
          init["frame"] = frame;
          init["a"] = a;
          init["e"] = e;
          init["i"] = inc;
          init["Omega"] = raan0 + p * 360.0 / n_planes;
          init["omega"] = omega;
          init["M"] = m0 + s * 360.0 / per_plane + p * phase;
          AddSatellite(Config(shared), init, fmt::format("{}-{}", prefix, idx++));
        }
      }
    } else {
      LUPNT_CHECK(false, "LunarNavConstellation needs a `satellites:` list or a `walker:` block",
                  "LunarNavConstellation");
    }
    Logger::Debug(fmt::format("LunarNavConstellation {} built {} satellites", name_, sats_.size()),
                  "LunarNavConstellation");
  }

  void LunarNavConstellation::Setup() {
    for (auto& sc : sats_) sc->SetSimulation(GetSimulation());
    Agent::Setup();  // schedule this constellation's own Step at frequency_ (if > 0)
  }

  void LunarNavConstellation::Step(Real t) {
    for (auto& sc : sats_) sc->Step(t);  // propagate + log each satellite's truth
    time_ = t;
  }

  Cart6 LunarNavConstellation::GetStateAt(Real t) const {
    if (sats_.empty()) return Cart6(Vec6::Zero(), Frame::MOON_CI);
    return sats_[0]->GetStateAt(t);
  }

  REGISTER_FACTORY_CLASS(Agent, LunarNavConstellation)

}  // namespace lupnt
