#pragma once

#include <string>
#include <vector>

#include "lupnt/applications/application.h"
#include "lupnt/core/config.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief Powered-descent guidance/control configuration: the truth descent trajectory of a
  /// lunar lander, in the local ENU tangent plane about the DEM site center. Either the built-in
  /// smoothstep descent (start/end East/North/Alt + heading) or an externally supplied reference
  /// trajectory (`ref_traj_enu`, e.g. from `pylupnt.lander_guidance`).
  struct LanderGncConfig {
    std::string start_epoch_utc = "2027-03-01T00:00:00";
    double duration_s = 300.0;  ///< total descent arc [s].
    double dt_s = 0.5;          ///< guidance step [s] (must match the nav app's dt_s).

    double descent_start_east_m = -2500.0;
    double descent_start_north_m = 400.0;
    double descent_end_east_m = 0.0;
    double descent_end_north_m = 0.0;
    double descent_start_alt_m = 2000.0;
    double descent_end_alt_m = 15.0;
    double descent_heading_deg = 0.0;

    /// Optional externally-supplied reference (truth) trajectory, `N x 3` rows of local
    /// East-North-Up position [m] about the DEM site center. When non-empty it overrides the
    /// built-in smoothstep descent and defines `N`.
    MatXd ref_traj_enu;
  };

  /// @brief Guidance/control application for a lunar lander, **hosted on a `Lander` agent**
  /// alongside a `LanderNavApp` (both in the agent's `applications:` list, guidance first).
  ///
  /// It owns the powered-descent **truth** trajectory: it builds the Moon-fixed position,
  /// velocity, attitude, and the body-frame specific force / angular rate from the ENU descent
  /// path (or a supplied reference trajectory), writes the host lander's truth state each epoch,
  /// and exposes the truth for the co-hosted navigation app to synthesize its measurements from.
  /// It reads only the shared `World` (terrain / ENU frame) — never the filter's estimate — so a
  /// closed-loop guidance law could later consume the nav estimate without changing this split.
  class LanderGncApp : public Application {
  public:
    LanderGncApp() = default;
    explicit LanderGncApp(Config& config);

    /// @brief Supply an externally-generated reference (truth) trajectory, `N x 3` ENU rows [m]
    /// about the DEM site center. Overrides the built-in smoothstep descent and sets `N`. Call
    /// before the `Simulation` runs.
    void SetReferenceTrajectoryEnu(const MatXd& ref_traj_enu) { cfg_.ref_traj_enu = ref_traj_enu; }

    void Setup() override;
    void Step(Real t) override;
    void Log(Real /*t*/) override {}

    /// @brief Build the truth trajectory if not already built (idempotent). Called on the first
    /// `Step`, and by a co-hosted `LanderNavApp` so it can read the truth regardless of app order.
    void EnsureInitialized();

    // ---- Truth accessors (valid after EnsureInitialized) ----------------------
    int N() const { return N_; }
    double dt() const { return dt_; }
    /// Moon-fixed (MOON_PA) truth position [m] at epoch k.
    const Vec3d& TruthPos(int k) const { return r_truth_[k]; }
    /// Moon-fixed (MOON_PA) truth velocity [m/s] at epoch k.
    const Vec3d& TruthVel(int k) const { return v_truth_[k]; }
    /// Body-to-MOON_PA truth attitude DCM at epoch k.
    const Mat3d& TruthAtt(int k) const { return R_truth_[k]; }
    /// Body-frame truth specific force [m/s^2] at epoch k (gravity removed).
    const Vec3d& SpecificForce(int k) const { return f_body_truth_[k]; }
    /// Body-frame truth angular rate [rad/s] at epoch k.
    const Vec3d& AngularRate(int k) const { return w_body_truth_[k]; }
    /// Truth height above terrain [m] at epoch k.
    double Altitude(int k) const { return Alt_[k]; }
    /// Truth ENU position (East, North, Up) [m] at epoch k.
    Vec3d EnuTruth(int k) const { return Vec3d(Ee_[k], Nn_[k], Uu_[k]); }

    const LanderGncConfig& config() const { return cfg_; }
    /// [N x 3] truth ENU trajectory (guidance path), for plotting.
    const MatXd& traj_enu_truth() const { return traj_enu_truth_; }
    /// [N] truth height above terrain, for plotting.
    const VecXd& alt_truth() const { return alt_truth_res_; }

  private:
    void InitScenario();
    void WriteLanderState(int k);

    LanderGncConfig cfg_;
    bool initialized_ = false;
    int N_ = 0;
    double dt_ = 0.5;

    // Local ENU tangent frame (from the shared World).
    Mat3d R_enu2pa_ = Mat3d::Identity();
    Vec3d up_hat_pa_ = Vec3d::UnitZ();

    // Truth trajectory.
    std::vector<double> Ee_, Nn_, Uu_, Alt_;
    std::vector<Vec3d> r_truth_, v_truth_;
    std::vector<Mat3d> R_truth_;
    std::vector<Vec3d> f_body_truth_, w_body_truth_;

    MatXd traj_enu_truth_;
    VecXd alt_truth_res_;
  };

}  // namespace lupnt
