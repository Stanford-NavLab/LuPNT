#include "lupnt/applications/ground_station/ground_station_tracking_app.h"

#include <limits>
#include <utility>

#include "lupnt/applications/ground_station/ground_station_manager_app.h"
#include "lupnt/lupnt.h"
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

    // Frame-invariant range / range-rate w.r.t. the station's inertial state.
    Vec3d dr = (xt.head(3) - st_world.head(3)).cast<double>();
    Vec3d dv = (xt.tail(3) - st_world.tail(3)).cast<double>();
    double rho = dr.norm();
    double rho_dot = dr.dot(dv) / rho;

    double range = use_range_ ? rho + SampleNormal(0.0, range_sigma_m_, &noise_rng_).val()
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
    manager_->AddMeasurement(m);
  }

  void GroundStationTrackingApp::Log(Real /*t*/) {
    DataLogger::Log(fmt::format("{}/num_measurements", name_), static_cast<double>(meas_t_.size()));
  }

  REGISTER_FACTORY_CLASS(Application, GroundStationTrackingApp)

}  // namespace lupnt
