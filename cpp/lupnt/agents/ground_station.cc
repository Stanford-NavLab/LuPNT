/**
 * @file groundstation.cc
 * @author Stanford NAV Lab
 * @brief  GroundStation Agent
 * @version 0.1
 * @date 2024-12-04
 *
 * @copyright Copyright (c) 2024
 *
 */
#include "lupnt/agents/ground_station.h"

#include <initializer_list>
#include <stdexcept>

#include "lupnt/core/asset_factory.h"
#include "lupnt/core/error.h"
#include "lupnt/dynamics/surface_dynamics.h"
#include "lupnt/interfaces/yaml.h"

namespace lupnt {
  namespace {
    Real ReadRequiredReal(const Config& config, std::initializer_list<const char*> keys) {
      for (const char* key : keys) {
        if (config[key]) return config[key].as<Real>();
      }
      throw std::runtime_error(
          fmt::format("GroundStation config missing required field {}", *keys.begin()));
    }

    Real ReadOptionalReal(const Config& config, std::initializer_list<const char*> keys,
                          Real fallback) {
      for (const char* key : keys) {
        if (config[key]) return config[key].as<Real>();
      }
      return fallback;
    }

    BodyId ReadBodyId(const Config& config, BodyId fallback) {
      if (config["body_id"]) return static_cast<BodyId>(config["body_id"].as<int>());
      if (!config["body"]) return fallback;

      std::string body = config["body"].as<std::string>();
      auto body_id = enum_cast<BodyId>(body);
      if (!body_id) throw std::runtime_error(fmt::format("Unknown body '{}'", body));
      return *body_id;
    }

    spice::GroundStationSpiceData MakeGroundStationData(const std::string& name,
                                                        double latitude_deg, double longitude_deg,
                                                        double altitude_m,
                                                        BodyId body_id = BodyId::EARTH) {
      BodyData body_data = GetBodyData(body_id);
      LatLonAlt lla(Vec3(latitude_deg, longitude_deg, altitude_m), body_data.fixed_frame);
      Cart3 position = LatLonAltToCart(lla, body_data.R, body_data.flattening);

      spice::GroundStationSpiceData data;
      data.name = name;
      data.body_id = body_id;
      data.frame = enum_name(body_data.fixed_frame);
      data.position_m = position.head(3).cast<double>();
      data.latitude_deg = latitude_deg;
      data.longitude_deg = longitude_deg;
      data.altitude_m = altitude_m;
      return data;
    }
  }  // namespace

  GroundStation::GroundStation(Config& config) : AgentWithDynamics() {
    ConfigureFromConfig(config);
  }

  GroundStation::GroundStation(const spice::GroundStationSpiceData& gs_data) : AgentWithDynamics() {
    ConfigureFromData(gs_data);
  }

  GroundStation GroundStation::FromSpice(Real t_tai, const std::string& station_name, BodyId center,
                                         const std::string& ref_frame,
                                         const std::string& ab_correction) {
    return GroundStation(
        spice::GetGroundStationDataSpice(t_tai, station_name, center, ref_frame, ab_correction));
  }

  GroundStation GroundStation::FromSpice(Real t_tai, int station_id, BodyId center,
                                         const std::string& ref_frame,
                                         const std::string& ab_correction) {
    return GroundStation(
        spice::GetGroundStationDataSpice(t_tai, station_id, center, ref_frame, ab_correction));
  }

  void GroundStation::ConfigureFromConfig(Config& config) {
    config_ = config;
    if (config["name"])
      SetName(config["name"].as<std::string>());
    else
      SetName(GetId());
    if (config["frequency"]) SetFrequency(config["frequency"].as<Real>());

    if (config["spice_kernel"]) {
      spice::LoadSpiceKernel(config["spice_kernel"].as<std::string>());
    }
    if (config["spice_kernels"]) {
      for (const auto& kernel : config["spice_kernels"]) {
        spice::LoadSpiceKernel(kernel.as<std::string>());
      }
    }

    BodyId body_id = ReadBodyId(config, BodyId::EARTH);
    Real epoch = ReadOptionalReal(config, {"epoch", "t_tai"}, 0.0);
    std::string ref_frame
        = config["spice_frame"] ? config["spice_frame"].as<std::string>() : "ITRF93";
    std::string ab_correction
        = config["ab_correction"] ? config["ab_correction"].as<std::string>() : "NONE";

    if (config["spice_name"]) {
      auto gs_data = spice::GetGroundStationDataSpice(epoch, config["spice_name"].as<std::string>(),
                                                      body_id, ref_frame, ab_correction);
      ConfigureFromData(gs_data);
      return;
    }

    if (config["spice_id"]) {
      auto gs_data = spice::GetGroundStationDataSpice(epoch, config["spice_id"].as<int>(), body_id,
                                                      ref_frame, ab_correction);
      ConfigureFromData(gs_data);
      return;
    }

    Real latitude_deg = ReadRequiredReal(config, {"latitude_deg", "latitude"});
    Real longitude_deg = ReadRequiredReal(config, {"longitude_deg", "longitude"});
    Real altitude_m = ReadRequiredReal(config, {"altitude_m", "altitude"});
    ConfigureFromData(MakeGroundStationData(GetName(), latitude_deg.val(), longitude_deg.val(),
                                            altitude_m.val(), body_id));

    // Optional onboard application (e.g. GroundStationOdtsApp). GroundStation builds
    // itself via ConfigureFromConfig rather than the Agent(Config&) constructor, so it
    // must opt into application creation explicitly.
    CreateApplication(config);
  }

  void GroundStation::ConfigureFromData(const spice::GroundStationSpiceData& gs_data) {
    ConfigureFromCartesian(gs_data.name, gs_data.position_m.cast<Real>(), gs_data.body_id);
    latitude_ = gs_data.latitude_deg;
    longitude_ = gs_data.longitude_deg;
    altitude_ = gs_data.altitude_m;
    naif_id_ = gs_data.naif_id;
    spice_frame_ = gs_data.frame;
  }

  void GroundStation::ConfigureFromCartesian(const std::string& name, const Vec3& position_m,
                                             BodyId body_id) {
    SetName(name);
    body_id_ = body_id;
    SetBodyId(body_id_);

    BodyData body_data = GetBodyData(body_id_);
    frame_ = body_data.fixed_frame;

    Cart3 position(position_m, frame_);
    LatLonAlt lla = CartToLatLonAlt(position, body_data.R, body_data.flattening);
    latitude_ = lla(0);
    longitude_ = lla(1);
    altitude_ = lla(2);

    Cart6 state(position_m, Vec3::Zero(), frame_);
    SetState(state);
    SetDynamics(MakePtr<StaticDynamics>());
    SetAttitude(Attitude(Vec4(1.0, 0.0, 0.0, 0.0), Vec3::Zero(), frame_));
  }

  REGISTER_FACTORY_CLASS(Agent, GroundStation)

}  // namespace lupnt
