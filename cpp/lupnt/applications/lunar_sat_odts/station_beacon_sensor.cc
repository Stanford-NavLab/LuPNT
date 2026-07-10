#include "lupnt/applications/lunar_sat_odts/station_beacon_sensor.h"

#include "lupnt/agents/spacecraft.h"
#include "lupnt/applications/lunar_sat_odts/ground_odts_app.h"
#include "lupnt/lupnt.h"
#include "lupnt/simulations/simulation.h"

namespace lupnt {

  StationBeaconSensor::StationBeaconSensor(Config& config) : Application(config) {
    auto g = [&](const char* k, double d) { return config[k].as<double>(d); };
    manager_name_ = config["manager"].as<std::string>();
    lat_deg_ = g("latitude_deg", 0.0);
    lon_deg_ = g("longitude_deg", 0.0);
    alt_m_ = g("altitude_m", 0.0);
    elevation_mask_deg_ = g("elevation_mask_deg", 5.0);
    pseudorange_sigma_m_ = g("pseudorange_sigma_m", 1.0);
    include_doppler_ = config["enable_station_doppler"].as<bool>(true);
    doppler_sigma_mps_ = g("station_doppler_sigma_mps", 1.0e-3);
    seed_ = config["seed"].as<int>(42);
    if (config["satellites"])
      for (const auto& s : config["satellites"]) sat_names_.push_back(s.as<std::string>());
  }

  void StationBeaconSensor::Setup() {
    LUPNT_CHECK(agent_, "Agent not set", "StationBeaconSensor");
    Application::Setup();  // schedule Step at APPLICATION priority (before the DEVICE-priority
                           // filter)
  }

  void StationBeaconSensor::Initialize() {
    Simulation* sim = agent_->GetSimulation();
    // Resolve the manager estimator.
    Agent* mgr = sim->GetAgent(manager_name_);
    LUPNT_CHECK(mgr, fmt::format("manager `{}` not found", manager_name_), "StationBeaconSensor");
    manager_ = dynamic_cast<GroundOdtsApp*>(mgr->GetApplication().get());
    LUPNT_CHECK(manager_, "StationBeaconSensor `manager` must host a GroundOdtsApp",
                "StationBeaconSensor");

    // Resolve tracked satellites and their filter-state block indices.
    sats_.clear();
    sat_index_.clear();
    for (const auto& sn : sat_names_) {
      auto* s = dynamic_cast<Spacecraft*>(sim->GetAgent(sn));
      LUPNT_CHECK(s, fmt::format("`{}` is not a Spacecraft", sn), "StationBeaconSensor");
      int idx = manager_->SatIndex(s->GetName());
      LUPNT_CHECK(idx >= 0, fmt::format("`{}` unknown to the manager", sn), "StationBeaconSensor");
      sats_.push_back(s);
      sat_index_.push_back(idx);
    }

    // Fixed station position in MOON_PA [m] from lat/lon/alt.
    VecX lla(3);
    lla << Real(lat_deg_), Real(lon_deg_), Real(alt_m_);
    State cart = LatLonAltToCart(lla, R_MOON, 0.0);
    r_bf_ = Vec3(cart.head(3));

    rng_.seed(static_cast<unsigned int>(seed_)
              + std::hash<std::string>{}(agent_->GetName()) % 99991u);
  }

  void StationBeaconSensor::Step(Real t) {
    if (!initialized_) {
      Initialize();
      initialized_ = true;
    }
    Real epoch_abs = GetLupntEpoch() + t;

    // Station MOON_CI position/velocity at this epoch.
    Vec6 st6;
    st6 << r_bf_, Vec3::Zero();
    Vec6 st_mci = ConvertFrame(epoch_abs, st6, Frame::MOON_PA, Frame::MOON_CI);
    Vec3d rst = Vec3(st_mci.head(3)).cast<double>();
    Vec3d vst = Vec3(st_mci.tail(3)).cast<double>();

    std::normal_distribution<double> nd(0.0, 1.0);
    for (size_t j = 0; j < sats_.size(); ++j) {
      VecXd xj = sats_[j]->GetTruthStateAt(t).cast<double>();
      // Elevation gate in the station's local MOON_PA frame.
      Vec6 xt;
      xt << xj.head(3).cast<Real>(), xj.segment(3, 3).cast<Real>();
      Vec6 xt_bf = ConvertFrame(epoch_abs, xt, Frame::MOON_CI, Frame::MOON_PA);
      Cart3 r_sat_bf(Vec3(xt_bf.head(3)), Frame::MOON_PA);
      Cart3 r_gs_bf(r_bf_, Frame::MOON_PA);
      State aer = CartToAzElRange(r_sat_bf, r_gs_bf);
      if ((aer(1) * DEG).val() <= elevation_mask_deg_) continue;

      Vec3d rj = xj.head(3), vj = xj.segment(3, 3);
      double bj = xj(6), dj = xj(7);
      IslStationObs m;
      m.t = t.val();
      m.sat_index = sat_index_[j];
      m.station_mci = rst;
      m.station_vel_mci = vst;
      m.pseudorange_m = (rj - rst).norm() + C * bj + pseudorange_sigma_m_ * nd(rng_);
      if (include_doppler_) {
        Vec3d u = (rj - rst) / (rj - rst).norm();
        m.has_doppler = true;
        m.doppler_mps = u.dot(vj - vst) + C * dj + doppler_sigma_mps_ * nd(rng_);
      }
      manager_->AddMeasurement(m);
    }
  }

  REGISTER_FACTORY_CLASS(Application, StationBeaconSensor)

}  // namespace lupnt
