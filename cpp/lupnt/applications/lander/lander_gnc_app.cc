#include "lupnt/applications/lander/lander_gnc_app.h"

#include <algorithm>
#include <cmath>

#include "lupnt/agents/agent.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/error.h"
#include "lupnt/simulations/simulation.h"
#include "lupnt/simulations/world.h"

namespace lupnt {

  namespace {
    // Smoothstep S(tau) = tau^2 (3 - 2 tau), clamped to [0, 1].
    double SmoothStep(double tau) {
      tau = std::clamp(tau, 0.0, 1.0);
      return tau * tau * (3.0 - 2.0 * tau);
    }
  }  // namespace

  LanderGncApp::LanderGncApp(Config& config) : Application(config) {
    LanderGncConfig c;
    c.start_epoch_utc = config["start_epoch_utc"].as<std::string>(c.start_epoch_utc);
    c.duration_s = config["duration_s"].as<double>(c.duration_s);
    c.dt_s = config["dt_s"].as<double>(c.dt_s);
    c.descent_start_east_m = config["descent_start_east_m"].as<double>(c.descent_start_east_m);
    c.descent_start_north_m = config["descent_start_north_m"].as<double>(c.descent_start_north_m);
    c.descent_end_east_m = config["descent_end_east_m"].as<double>(c.descent_end_east_m);
    c.descent_end_north_m = config["descent_end_north_m"].as<double>(c.descent_end_north_m);
    c.descent_start_alt_m = config["descent_start_alt_m"].as<double>(c.descent_start_alt_m);
    c.descent_end_alt_m = config["descent_end_alt_m"].as<double>(c.descent_end_alt_m);
    c.descent_heading_deg = config["descent_heading_deg"].as<double>(c.descent_heading_deg);
    cfg_ = c;
  }

  void LanderGncApp::Setup() {
    LUPNT_CHECK(agent_, "Agent not set", "LanderGncApp::Setup");
    dt_ = cfg_.dt_s;
    // Defer the trajectory build to the first Step so a reference trajectory set via
    // SetReferenceTrajectoryEnu (after construction, before run()) is honored.
    Simulation* sim = agent_->GetSimulation();
    sim->Schedule(0.0, [this](Real t) { Step(t); }, 1.0 / dt_, Event::Priority::APPLICATION);
  }

  void LanderGncApp::EnsureInitialized() {
    if (initialized_) return;
    InitScenario();
    initialized_ = true;
  }

  void LanderGncApp::InitScenario() {
    World* world = agent_->GetWorld();
    LUPNT_CHECK(world && world->HasTerrain(),
                "LanderGncApp requires a World with a `dem:` terrain block", "LanderGncApp::Setup");

    dt_ = cfg_.dt_s;
    const double dt = dt_;

    const bool use_ref_traj = cfg_.ref_traj_enu.rows() > 0;
    LUPNT_CHECK(!use_ref_traj || cfg_.ref_traj_enu.cols() == 3,
                "ref_traj_enu must have 3 columns (East, North, Up)", "LanderGncApp::Setup");
    N_ = use_ref_traj ? static_cast<int>(cfg_.ref_traj_enu.rows())
                      : std::max(2, static_cast<int>(std::round(cfg_.duration_s / cfg_.dt_s)) + 1);

    // ---- Local ENU terrain frame (from the shared World) ----------------------
    R_enu2pa_ = world->REnuToWorld();
    up_hat_pa_ = R_enu2pa_.col(2);

    auto EnuToPa
        = [&](double E, double Nn, double U) -> Vec3d { return world->EnuToWorld(E, Nn, U); };
    auto Elevation = [&](double E, double Nn) -> double { return world->GetElevation(E, Nn); };
    auto BodyToPa = [&](double hdg) -> Mat3d {
      Vec3d xb(std::cos(hdg), std::sin(hdg), 0.0);
      Vec3d zb(0.0, 0.0, 1.0);
      Vec3d yb = zb.cross(xb);
      Mat3d R_body2enu;
      R_body2enu.col(0) = xb;
      R_body2enu.col(1) = yb;
      R_body2enu.col(2) = zb;
      return R_enu2pa_ * R_body2enu;
    };

    // ---- Descent truth: ENU path, Moon-fixed pos/vel/accel, attitude/rate ------
    const double hdg = cfg_.descent_heading_deg * RAD;
    Ee_.assign(N_, 0.0);
    Nn_.assign(N_, 0.0);
    Uu_.assign(N_, 0.0);
    Alt_.assign(N_, 0.0);
    r_truth_.assign(N_, Vec3d::Zero());
    v_truth_.assign(N_, Vec3d::Zero());
    R_truth_.assign(N_, Mat3d::Identity());
    for (int k = 0; k < N_; ++k) {
      if (use_ref_traj) {
        Ee_[k] = cfg_.ref_traj_enu(k, 0);
        Nn_[k] = cfg_.ref_traj_enu(k, 1);
        Uu_[k] = cfg_.ref_traj_enu(k, 2);
        Alt_[k] = Uu_[k] - Elevation(Ee_[k], Nn_[k]);
      } else {
        double tau = (cfg_.duration_s > 0.0) ? (k * dt / cfg_.duration_s) : 1.0;
        double s = SmoothStep(tau);
        Ee_[k]
            = cfg_.descent_start_east_m + (cfg_.descent_end_east_m - cfg_.descent_start_east_m) * s;
        Nn_[k] = cfg_.descent_start_north_m
                 + (cfg_.descent_end_north_m - cfg_.descent_start_north_m) * s;
        Alt_[k]
            = cfg_.descent_start_alt_m + (cfg_.descent_end_alt_m - cfg_.descent_start_alt_m) * s;
        Uu_[k] = Elevation(Ee_[k], Nn_[k]) + Alt_[k];
      }
      R_truth_[k] = BodyToPa(hdg);
    }
    for (int k = 0; k < N_; ++k) r_truth_[k] = EnuToPa(Ee_[k], Nn_[k], Uu_[k]);
    for (int k = 0; k < N_; ++k) {
      int kp = std::min(k + 1, N_ - 1), km = std::max(k - 1, 0);
      double span = (kp - km) * dt;
      v_truth_[k]
          = (span > 0.0) ? Vec3d((r_truth_[kp] - r_truth_[km]) / span) : Vec3d(Vec3d::Zero());
    }
    f_body_truth_.assign(N_, Vec3d::Zero());
    w_body_truth_.assign(N_, Vec3d::Zero());
    for (int k = 0; k < N_; ++k) {
      int kp = std::min(k + 1, N_ - 1), km = std::max(k - 1, 0);
      double span = (kp - km) * dt;
      Vec3d a_total
          = (span > 0.0) ? Vec3d((v_truth_[kp] - v_truth_[km]) / span) : Vec3d(Vec3d::Zero());
      f_body_truth_[k] = R_truth_[k].transpose() * (a_total - world->Gravity(r_truth_[k]));
      Mat3d Rdot
          = (span > 0.0) ? Mat3d((R_truth_[kp] - R_truth_[km]) / span) : Mat3d(Mat3d::Zero());
      Mat3d Wx = R_truth_[k].transpose() * Rdot;
      w_body_truth_[k] = Vec3d(Wx(2, 1), Wx(0, 2), Wx(1, 0));
    }

    // Truth-trajectory result series (for plotting the guidance path).
    traj_enu_truth_ = MatXd::Zero(N_, 3);
    alt_truth_res_ = VecXd::Zero(N_);
    for (int k = 0; k < N_; ++k) {
      traj_enu_truth_.row(k) = Vec3d(Ee_[k], Nn_[k], Uu_[k]).transpose();
      alt_truth_res_(k) = Alt_[k];
    }

    // Seed the host lander's truth state at epoch 0.
    WriteLanderState(0);
  }

  void LanderGncApp::WriteLanderState(int k) {
    auto* lander = dynamic_cast<AgentWithDynamics*>(agent_);
    if (!lander) return;
    Vec6 rv;
    rv << r_truth_[k].cast<Real>(), v_truth_[k].cast<Real>();
    lander->SetTime(k * dt_);
    lander->SetState(Cart6(rv, Frame::MOON_PA));
  }

  void LanderGncApp::Step(Real t) {
    EnsureInitialized();
    int k = static_cast<int>(std::lround(t.val() / dt_));
    if (k < 1 || k >= N_) return;
    WriteLanderState(k);
  }

  REGISTER_FACTORY_CLASS(Application, LanderGncApp)

}  // namespace lupnt
