#include "lupnt/applications/isl_odts_app.h"

#include <cmath>

#include "lupnt/core/error.h"
#include "lupnt/lupnt.h"
#include "lupnt/measurements/measurement.h"
#include "lupnt/measurements/measurement_utils.h"
#include "lupnt/numerics/filters/filter_utils.h"

namespace lupnt {

  namespace {

    constexpr int kSubStateSize = IslOdtsApp::kSubStateSize;

    // Extract sub-state `x.segment(offset, 8)` as a properly labeled
    // JointOrbitClockState, suitable for propagation by a JointOrbitClockDynamics.
    JointOrbitClockState ExtractSubState(const State& x, int offset) {
      State sub(kSubStateSize);
      sub = x.segment(offset, kSubStateSize);
      sub.SetFrame(x.GetFrame());
      return JointOrbitClockState(sub);
    }

    // Block-diagonal [own(8), consider_1(8), ..., consider_L(8)] dynamics: every
    // block is propagated independently (no dynamical coupling between satellites),
    // by the same (deterministic, filter-fidelity) JointOrbitClockDynamics model.
    FilterDynamicsFunction MakeCrosslinkDynamicsFunction(Ptr<JointOrbitClockDynamics> dyn,
                                                         int n_sat) {
      return [dyn, n_sat](const State& x, Real t0, Real tf, const State*, MatXd* F) -> State {
        const int n = n_sat * kSubStateSize;
        State xf(n);
        xf.SetFrame(x.GetFrame());
        if (F != nullptr) F->setZero(n, n);
        for (int j = 0; j < n_sat; ++j) {
          const int off = j * kSubStateSize;
          JointOrbitClockState xj = ExtractSubState(x, off);
          if (F != nullptr) {
            MatXd Fj;
            State xfj = dyn->Propagate(xj, t0, tf, nullptr, &Fj);
            xf.segment(off, kSubStateSize) = xfj;
            F->block(off, off, kSubStateSize, kSubStateSize) = Fj;
          } else {
            State xfj = dyn->Propagate(xj, t0, tf, nullptr);
            xf.segment(off, kSubStateSize) = xfj;
          }
        }
        return xf;
      };
    }

    MatXd ProcessNoiseSubState(double sigma_a_mps2, Real dt) {
      MatXd Q = MatXd::Zero(kSubStateSize, kSubStateSize);
      Mat3d Q_acc = std::pow(sigma_a_mps2, 2) * Mat3d::Identity();
      Q.block(0, 0, 6, 6) = ProcessNoisePosVel(Q_acc, dt);
      Q.block(6, 6, 2, 2) = ClockDynamics::TwoStateNoise(ClockModel::OCXO, dt).cast<double>();
      return Q;
    }

    ProcessNoiseFunction MakeCrosslinkProcessNoiseFunction(double sigma_a_mps2, int n_sat) {
      return [sigma_a_mps2, n_sat](const State& x, Real t0, Real tf) {
        const double dt = std::abs((tf - t0).val());
        MatXd Q_sub = ProcessNoiseSubState(sigma_a_mps2, dt);
        MatXd Q = MatXd::Zero(x.size(), x.size());
        for (int j = 0; j < n_sat; ++j) {
          Q.block(j * kSubStateSize, j * kSubStateSize, kSubStateSize, kSubStateSize) = Q_sub;
        }
        return Q;
      };
    }

  }  // namespace

  IslOdtsApp::IslOdtsApp(IslOdtsAppParams params)
      : LunaNetSubApp("isl_odts"), params_(std::move(params)) {}

  void IslOdtsApp::Configure(Real t0, const VecXd& x0, const MatXd& P0,
                             Ptr<JointOrbitClockDynamics> filter_dynamics) {
    t0_ = t0;
    x0_ = x0;
    P0_ = P0;
    filter_dynamics_ = std::move(filter_dynamics);
  }

  void IslOdtsApp::Setup(LunaNetSatApp& app) {
    LunaNetSubApp::Setup(app);
    LUPNT_CHECK(filter_dynamics_ != nullptr, "Call Configure() before Setup()", "IslOdtsApp");
    LUPNT_CHECK(params_.n_sat >= 2, "IslOdtsAppParams.n_sat must be >= 2", "IslOdtsApp");

    const int n_sat = params_.n_sat;
    const int n_links = n_sat - 1;
    filter_ = MakePtr<SchmidtEKF>(kSubStateSize * n_links);  // trailing consider blocks

    State x0_state(kSubStateSize * n_sat);
    x0_state = x0_.cast<Real>();
    x0_state.SetFrame(Frame::MOON_CI);
    filter_->SetTime(t0_);
    filter_->SetState(x0_state);
    filter_->SetCovariance(P0_);
    filter_->SetDynamicsFunction(MakeCrosslinkDynamicsFunction(filter_dynamics_, n_sat));
    filter_->SetProcessNoiseFunction(
        MakeCrosslinkProcessNoiseFunction(params_.process_accel_sigma_mps2, n_sat));
    // The measurement function is (re)installed every epoch in Step() because its
    // tracked-GPS set changes.
    filter_->SetOutlierThreshold(params_.outlier_threshold);
  }

  void IslOdtsApp::StageMeasurements(const IslOdtsMeasurementEpoch& meas) {
    staged_ = meas;
    has_staged_ = true;
  }

  void IslOdtsApp::Step(Real t) {
    LUPNT_CHECK(filter_ != nullptr, "Call Setup() before Step()", "IslOdtsApp");
    LUPNT_CHECK(has_staged_, "StageMeasurements() must be called before each Step()", "IslOdtsApp");
    const int n_links = params_.n_sat - 1;
    const int n_anchor = static_cast<int>(staged_.anchor_pos_mci.size());

    filter_->Predict(t);

    // Point the hub filter at this epoch's measurement geometry: the crosslinks plus
    // (when an anchor is staged) any one-way pseudoranges to known-position anchors
    // (e.g. a lunar surface station). With no staged anchor this is the crosslink-only
    // model. The model itself lives in `IslCrosslinkMeasurement`; here we only supply
    // this epoch's config/geometry.
    IslCrosslinkMeasurement::Config meas_cfg;
    meas_cfg.sub_state_size = kSubStateSize;
    meas_cfg.n_links = n_links;
    meas_cfg.anchor_pos_mci = staged_.anchor_pos_mci;
    meas_cfg.sigma_range_m = params_.range_sigma_m;
    meas_cfg.sigma_range_rate_mps = params_.range_rate_sigma_mps;
    meas_cfg.sigma_pseudorange_m = params_.pseudorange_sigma_m;
    filter_->SetMeasurementFunction(IslCrosslinkMeasurement(meas_cfg).CreateFunction());

    VecXd y_obs(2 * n_links + n_anchor);
    for (int i = 0; i < n_links; ++i) {
      y_obs(2 * i) = staged_.crosslink_range_m(i);
      y_obs(2 * i + 1) = staged_.crosslink_range_rate_mps(i);
    }
    for (int a = 0; a < n_anchor; ++a) y_obs(2 * n_links + a) = staged_.anchor_pseudorange_m(a);

    filter_->Update(y_obs);
    prefit_resid_ = filter_->GetMeasurementResidual();
    has_staged_ = false;
  }

  void IslOdtsApp::Finish() { LunaNetSubApp::Finish(); }

}  // namespace lupnt
