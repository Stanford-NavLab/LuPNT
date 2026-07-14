#include "lupnt/applications/ground_station/ground_station_tracking_app.h"

#include <limits>
#include <utility>

#include "lupnt/agents/ground_station.h"
#include "lupnt/applications/ground_station/ground_station_manager_app.h"
#include "lupnt/lupnt.h"
#include "lupnt/measurements/ground_station_corrections.h"
#include "lupnt/simulations/world.h"

namespace lupnt {

  GroundStationTrackingApp::GroundStationTrackingApp(Config& config) : Application(config) {
    Logger::Debug(fmt::format("Creating {}", name_), "GroundStationTrackingApp");

    LUPNT_CHECK(config["target"], "GroundStationTrackingApp requires a `target` agent name",
                "GroundStationTrackingApp");
    target_name_ = config["target"].as<std::string>();
    LUPNT_CHECK(config["manager"], "GroundStationTrackingApp requires a `manager` agent name",
                "GroundStationTrackingApp");
    manager_name_ = config["manager"].as<std::string>();

    elevation_mask_deg_ = config["elevation_mask_deg"].as<double>(elevation_mask_deg_);
    use_range_ = config["use_range"].as<bool>(use_range_);
    use_range_rate_ = config["use_range_rate"].as<bool>(use_range_rate_);
    range_sigma_m_ = config["range_sigma_m"].as<double>(range_sigma_m_);
    range_rate_sigma_mps_ = config["range_rate_sigma_mps"].as<double>(range_rate_sigma_mps_);
    seed_ = config["seed"].as<int>(seed_);

    // Optional truth-model corrections (Earth signal-path delays + station-location tide).
    apply_solid_earth_tide_ = config["apply_solid_earth_tide"].as<bool>(apply_solid_earth_tide_);
    apply_troposphere_ = config["apply_troposphere"].as<bool>(apply_troposphere_);
    apply_ionosphere_ = config["apply_ionosphere"].as<bool>(apply_ionosphere_);
    apply_shapiro_ = config["apply_shapiro"].as<bool>(apply_shapiro_);
    tropo_pressure_hpa_ = config["troposphere_pressure_hpa"].as<double>(tropo_pressure_hpa_);
    tropo_temperature_k_ = config["troposphere_temperature_k"].as<double>(tropo_temperature_k_);
    tropo_humidity_pct_ = config["troposphere_humidity_pct"].as<double>(tropo_humidity_pct_);
    iono_vtec_tecu_ = config["ionosphere_vtec_tecu"].as<double>(iono_vtec_tecu_);
    signal_frequency_hz_ = config["signal_frequency_hz"].as<double>(signal_frequency_hz_);
    corrections_enabled_
        = apply_solid_earth_tide_ || apply_troposphere_ || apply_ionosphere_ || apply_shapiro_;

    LUPNT_CHECK(use_range_ || use_range_rate_,
                "GroundStationTrackingApp needs at least one of use_range / use_range_rate",
                "GroundStationTrackingApp");
  }

  void GroundStationTrackingApp::Setup() {
    LUPNT_CHECK(agent_, "Agent not set", "GroundStationTrackingApp");
    Simulation* sim = agent_->GetSimulation();
    epoch0_ = GetLupntEpoch();

    // This station's fixed position in its body-fixed frame (velocity zero there).
    auto* gs = dynamic_cast<AgentWithDynamics*>(agent_);
    LUPNT_CHECK(gs, "GroundStationTrackingApp must run on an AgentWithDynamics (GroundStation)",
                "GroundStationTrackingApp");
    station_name_ = gs->GetName();
    State gs_state = gs->GetState();
    station_r_ = gs_state.head(3);
    station_frame_ = gs_state.GetFrame();

    // Geodetic latitude/height for the troposphere and solid-tide models (if requested).
    if (auto* gs_geo = dynamic_cast<GroundStation*>(agent_)) {
      station_lat_rad_ = gs_geo->GetLatitudeDegDouble() * RAD;
      station_height_m_ = gs_geo->GetAltitudeMDouble();
    }

    // Resolve the manager app and register this station with it.
    Agent* mgr_agent = sim->GetAgent(manager_name_);
    manager_ = dynamic_cast<GroundStationManagerApp*>(mgr_agent->GetApplication().get());
    LUPNT_CHECK(manager_,
                fmt::format("Manager agent `{}` has no GroundStationManagerApp", manager_name_),
                "GroundStationTrackingApp");
    station_id_ = manager_->RegisterStation(station_name_);

    noise_rng_.seed(static_cast<unsigned int>(seed_) + static_cast<unsigned int>(station_id_));

    Logger::Info(fmt::format("{}: tracking `{}` -> manager `{}` (elevation mask {:.1f} deg)", name_,
                             target_name_, manager_name_, elevation_mask_deg_),
                 "GroundStationTrackingApp");

    // Schedule periodic measurement Steps (base class, using frequency_).
    Application::Setup();
  }

  void GroundStationTrackingApp::Step(Real t) {
    World* world = agent_->GetWorld();
    LUPNT_CHECK(world, "GroundStationTrackingApp requires a World (define a `world:` block)",
                "GroundStationTrackingApp");
    Frame world_frame = world->GetFrame();
    Real epoch_abs = epoch0_ + t;

    // Target truth state and this station's state, both in the world (inertial) frame.
    Vec6 xt = world->GetStateAt(target_name_, t);
    Vec6 st6;
    st6 << station_r_, Vec3::Zero();
    Vec6 st_world = ConvertFrame(epoch_abs, st6, station_frame_, world_frame);

    // Topocentric elevation gate: convert the satellite to the station body-fixed frame.
    Vec6 xt_bf = ConvertFrame(epoch_abs, xt, world_frame, station_frame_);
    Cart3 r_sat_bf(Vec3(xt_bf.head(3)), station_frame_);
    Cart3 r_gs_bf(station_r_, station_frame_);
    State aer = CartToAzElRange(r_sat_bf, r_gs_bf);
    double elevation_deg = (aer(1) * DEG).val();

    elev_t_.push_back(t.val());
    elev_deg_.push_back(elevation_deg);
    if (elevation_deg <= elevation_mask_deg_) return;

    // Optional signal-path / station-location corrections. These act on the *truth*
    // geometry only: the nominal station position `st_world` is still reported to the
    // manager, so an estimator that models a pure geometric range sees the corrections as
    // realistic measurement errors. `st_truth` carries the solid-tide-displaced station.
    Vec6 st_truth = st_world;
    Vec3d tide_disp = Vec3d::Zero();  // solid-tide station displacement, world frame [m]
    double tropo_delay = 0.0, iono_delay = 0.0, shapiro_delay = 0.0;  // per-component [m]
    if (corrections_enabled_) {
      Vec3d earth_w = GetBodyPosVel(epoch_abs, BodyId::EARTH, world_frame).head(3).cast<double>();
      Vec3d sun_w = GetBodyPosVel(epoch_abs, BodyId::SUN, world_frame).head(3).cast<double>();

      if (apply_solid_earth_tide_) {
        // Geocentric station / tide-body vectors (world frame is Moon-centered, so the Moon
        // sits at the origin). The displacement is returned in the same inertial frame.
        Vec3d station_geo = st_world.head(3).cast<double>() - earth_w;
        Vec3d moon_geo = -earth_w;
        Vec3d sun_geo = sun_w - earth_w;
        tide_disp = SolidEarthTideDisplacement(
            station_geo, {{moon_geo, GM_MOON}, {sun_geo, GM_SUN}}, GM_EARTH, R_EARTH);
        st_truth.head(3) += tide_disp.cast<Real>();
      }
      if (apply_troposphere_) {
        double p_hpa = tropo_pressure_hpa_ > 0.0 ? tropo_pressure_hpa_
                                                 : StandardAtmospherePressureHPa(station_height_m_);
        double t_k = tropo_temperature_k_ > 0.0 ? tropo_temperature_k_
                                                : StandardAtmosphereTemperatureK(station_height_m_);
        tropo_delay
            = TroposphereDelaySaastamoinen(elevation_deg * RAD, station_lat_rad_, station_height_m_,
                                           p_hpa, t_k, tropo_humidity_pct_);
      }
      if (apply_ionosphere_) {
        iono_delay = IonosphereDelayThinShell(elevation_deg * RAD, iono_vtec_tecu_,
                                              signal_frequency_hz_, station_height_m_);
      }
      if (apply_shapiro_) {
        Vec3d r_tx = st_truth.head(3).cast<double>();
        Vec3d r_rx = xt.head(3).cast<double>();
        shapiro_delay = ShapiroRangeDelay(r_tx, r_rx, {{earth_w, GM_EARTH}, {sun_w, GM_SUN}});
      }
    }
    double path_delay = tropo_delay + iono_delay + shapiro_delay;  // total added to range [m]

    // Frame-invariant range / range-rate w.r.t. the (possibly tide-displaced) station state.
    Vec3d dr = (xt.head(3) - st_truth.head(3)).cast<double>();
    Vec3d dv = (xt.tail(3) - st_truth.tail(3)).cast<double>();
    double rho = dr.norm();
    double rho_dot = dr.dot(dv) / rho;

    double range = use_range_
                       ? rho + path_delay + SampleNormal(0.0, range_sigma_m_, &noise_rng_).val()
                       : std::numeric_limits<double>::quiet_NaN();
    double range_rate = use_range_rate_
                            ? rho_dot + SampleNormal(0.0, range_rate_sigma_mps_, &noise_rng_).val()
                            : std::numeric_limits<double>::quiet_NaN();

    meas_t_.push_back(t.val());
    meas_range_.push_back(range);
    meas_range_rate_.push_back(range_rate);

    // Report the observation to the centralized manager/estimator.
    StationMeasurement m;
    m.t = t.val();
    m.station_id = station_id_;
    m.epoch_index = manager_->EpochIndex(t.val());
    m.station_mci = st_world.cast<double>();
    m.has_range = use_range_;
    m.has_range_rate = use_range_rate_;
    m.range = range;
    m.range_rate = range_rate;
    m.range_sigma = range_sigma_m_;
    m.range_rate_sigma = range_rate_sigma_mps_;
    m.tropo_delay_m = tropo_delay;
    m.iono_delay_m = iono_delay;
    m.shapiro_delay_m = shapiro_delay;
    m.tide_disp_m = tide_disp;
    manager_->AddMeasurement(m);
  }

  void GroundStationTrackingApp::Log(Real /*t*/) {
    DataLogger::Log(fmt::format("{}/num_measurements", name_), static_cast<double>(meas_t_.size()));
  }

  REGISTER_FACTORY_CLASS(Application, GroundStationTrackingApp)

}  // namespace lupnt
