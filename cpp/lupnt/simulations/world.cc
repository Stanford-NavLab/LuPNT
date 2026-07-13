#include "lupnt/simulations/world.h"

#include "lupnt/agents/agent.h"
#include "lupnt/dynamics/numerical_orbit_dynamics.h"
#include "lupnt/lupnt.h"
#include "lupnt/simulations/simulation.h"

namespace lupnt {

  namespace {
    // Point-mass GM for a named central body (surface INS gravity). Defaults to the Moon.
    double GmForBody(const std::string& body) {
      if (body == "EARTH") return GM_EARTH;
      if (body == "SUN") return GM_SUN;
      if (body == "MOON") return GM_MOON;
      Logger::Warn(fmt::format("World: unknown gravity body `{}`, using MOON", body), "World");
      return GM_MOON;
    }
  }  // namespace

  World::World(Config& world_config) {
    // Epoch: the global LuPNT epoch, set by Simulation from the scenario `epoch:`.
    epoch_ = GetLupntEpoch();

    if (world_config["frame"]) {
      frame_ = enum_cast<Frame>(world_config["frame"].as<std::string>()).value();
    }

    // Orbital force model (optional).
    if (world_config["force_model"]) {
      has_force_model_ = true;
      force_model_ = world_config["force_model"];
    }

    // Ionosphere/plasmasphere signal-delay environment (optional; a shared truth
    // property read by the GNSS ODTS app, hence under `world:` not the receiver app).
    if (world_config["plasma"]) {
      has_plasma_ = true;
      plasma_ = world_config["plasma"];
    }

    // Central-body point-mass gravity for surface INS (optional; defaults to Moon).
    gm_ = GM_MOON;
    if (world_config["gravity"] && world_config["gravity"]["body"]) {
      gm_ = GmForBody(world_config["gravity"]["body"].as<std::string>());
    }

    // Terrain / DEM (optional, surface scenarios). The ENU tangent frame is defined
    // relative to the Moon body-fixed frame (MOON_PA), which is the world frame for
    // surface-navigation scenarios.
    if (world_config["dem"]) {
      const auto& d = world_config["dem"];
      site_lat_deg_ = d["site_lat_deg"].as<double>();
      site_lon_deg_ = d["site_lon_deg"].as<double>();
      double half_width_m = d["half_width_m"].as<double>(5000.0);
      double max_res_m = d["max_res_m"].as<double>(20.0);

      // Optional explicit tile: skip the PGDA download and crop this GeoTIFF directly. A
      // relative path is resolved against the current working directory. Used to run surface
      // scenarios (and their tests) against a small bundled fixture instead of the 41 MB tile.
      std::filesystem::path dem_file;
      if (d["dem_file"]) {
        dem_file = d["dem_file"].as<std::string>();
        if (!dem_file.empty() && dem_file.is_relative()) {
          std::error_code ec;
          dem_file = std::filesystem::absolute(dem_file, ec);
        }
      }

      dem_ = LoadLolaDem(site_lat_deg_, site_lon_deg_, half_width_m, max_res_m, dem_file);
      dem_cx_ = dem_.center_x();
      dem_cy_ = dem_.center_y();

      Vec3 r_center = LatLonAltToCart(Vec3(site_lat_deg_, site_lon_deg_, 0.0), R_MOON);
      Cart3 r_center_state(r_center, Frame::MOON_PA);
      R_enu2world_ = RotEastNorthUpToCart(r_center_state, R_MOON).cast<double>();
      r_center_world_ = r_center.cast<double>();
      has_dem_ = true;
      Logger::Debug(fmt::format("World: loaded DEM site {} ({}, {})", dem_.site().id, site_lat_deg_,
                                site_lon_deg_),
                    "World");
    }

    Logger::Debug(fmt::format("World created (epoch={} s past J2000)", epoch_.val()), "World");
  }

  Ptr<NBodyDynamics> World::MakeDynamics() const {
    LUPNT_CHECK(has_force_model_, "World has no `force_model:` block; MakeDynamics unavailable",
                "World");
    Config cfg = force_model_;  // NBodyDynamics(Config&) parses bodies/frame/SRP/integrator
    auto dynamics = MakePtr<NBodyDynamics>(cfg);
    dynamics->SetFrame(frame_);
    dynamics->SetAutodiff(true);  // enable the analytic state-transition matrix
    return dynamics;
  }

  Ptr<NBodyDynamics> World::MakeDynamics(double cr, double area_m2, double mass_kg) const {
    auto dynamics = MakeDynamics();
    dynamics->SetSrpCoefficient(cr, area_m2, mass_kg);
    return dynamics;
  }

  Ptr<NBodyDynamics> World::MakeTruthDynamics() const {
    auto dynamics = MakeDynamics();
    dynamics->SetAutodiff(false);  // truth propagation needs no state-transition matrix
    return dynamics;
  }

  Vec3d World::Gravity(const Vec3d& r) const {
    double rn = r.norm();
    return (rn > 0.0) ? Vec3d(-gm_ / (rn * rn * rn) * r) : Vec3d(Vec3d::Zero());
  }

  double World::GetElevation(double east_m, double north_m) const {
    LUPNT_CHECK(has_dem_, "World has no `dem:` terrain block", "World");
    return dem_.GetElevation(dem_cx_ + east_m, dem_cy_ + north_m);
  }

  Vec3d World::EnuToWorld(double east_m, double north_m, double up_m) const {
    LUPNT_CHECK(has_dem_, "World has no `dem:` terrain block", "World");
    return r_center_world_ + R_enu2world_ * Vec3d(east_m, north_m, up_m);
  }

  Cart6 World::GetStateAt(const std::string& agent_name, Real t) const {
    LUPNT_CHECK(sim_, "World has no Simulation registered (call SetWorld on the Simulation)",
                "World");
    Agent* agent = sim_->GetAgent(agent_name);
    return agent->GetStateAt(t);
  }

}  // namespace lupnt
