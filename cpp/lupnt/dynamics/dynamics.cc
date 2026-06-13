#include "lupnt/dynamics/dynamics.h"

#include "lupnt/core/asset_factory.h"
#include "lupnt/core/logger.h"
#include "lupnt/core/progress_bar.h"

namespace lupnt {

  Dynamics::Dynamics(Config& config) {
    if (config["print_progress"]) SetPrintProgress(config["print_progress"].as<bool>());
  }

  State Dynamics::Propagate(const State& x0, Real t0, Real tf, const State* u, MatXd* stm) {
    if (stm == nullptr) return Propagate(x0, t0, tf, u);
    auto func = [=, this](const VecX& x) {
      State xs = State(x, x0.GetNames(), x0.GetUnits());
      return Propagate(xs, t0, tf, u);
    };
    VecX xf_tmp;
    VecX x0_tmp = x0.cast<double>();
    jacobian(func, wrt(x0_tmp), at(x0_tmp), xf_tmp, *stm);
    State xf(xf_tmp, x0.GetNames(), x0.GetUnits());
    return xf;
  }

  MatX Dynamics::Propagate(const State& x0, const VecX& tfs, const State* u) {
    MatX xf = MatX::Zero(tfs.size(), x0.size());
    Ptr<ProgressBar> pbar = nullptr;
    if (print_progress_) pbar = Logger::GetProgressBar(tfs.size(), "Propagating", "Dynamics");
    xf.row(0) = x0;
    for (int i = 1; i < tfs.size(); i++) {
      Real t0_i = tfs(i - 1);
      Real tf_i = tfs(i);
      VecX x0_i = xf.row(i - 1);
      VecX xf_i = Propagate(x0_i, t0_i, tf_i, u);
      xf.row(i) = xf_i;
      if (print_progress_) pbar->Update(i);
    }
    if (print_progress_) pbar->Finish();
    return xf;
  }

  State Dynamics::PropagateWithParams(const State& x0, Real t0, Real tf, const ParamState& params,
                                      const State* u) {
    // Set parameters
    SetParams(params);
    return Propagate(x0, t0, tf, u);
  }

  State Dynamics::PropagateWithParams(const State& x0, Real t0, Real tf, const ParamState& params,
                                      const State* u, MatXd* stm_state, MatXd* stm_param) {
    if (stm_state == nullptr && stm_param == nullptr)
      return PropagateWithParams(x0, t0, tf, params, u);
    auto func = [=, this](const VecX& x, const VecX& p) {
      State xs = State(x, x0.GetNames(), x0.GetUnits());
      ParamState ps = ParamState(p, params.GetNames(), params.GetUnits());
      return PropagateWithParams(xs, t0, tf, ps, u);
    };
    VecX xf_tmp;
    VecX x0_tmp = x0.cast<double>();
    VecX p0_tmp = params.cast<double>();
    jacobian(func, wrt(x0_tmp), at(x0_tmp, p0_tmp), xf_tmp, *stm_state);
    jacobian(func, wrt(p0_tmp), at(x0_tmp, p0_tmp), xf_tmp, *stm_param);

    State xf(xf_tmp, x0.GetNames(), x0.GetUnits());

    return xf;
  }

  MatX Dynamics::PropagateWithParams(const State& x0, const VecX& tfs, const ParamState& params,
                                     const State* u) {
    // Set parameters
    MatX xf = MatX::Zero(tfs.size(), x0.size());
    Ptr<ProgressBar> pbar = nullptr;
    if (print_progress_) pbar = Logger::GetProgressBar(tfs.size(), "Propagating", "Dynamics");
    xf.row(0) = x0;
    for (int i = 1; i < tfs.size(); i++) {
      Real t0_i = tfs(i - 1);
      Real tf_i = tfs(i);
      VecX x0_i = xf.row(i - 1);
      VecX xf_i = PropagateWithParams(x0_i, t0_i, tf_i, params, u);
      xf.row(i) = xf_i;
      if (print_progress_) pbar->Update(i);
    }
    if (print_progress_) pbar->Finish();
    return xf;
  }

  // Define the GetRegistry function for this specialization (must come before explicit
  // instantiation)
  template <> std::unordered_map<std::string, AssetFactory<Dynamics, Config&>::Creator>&
  AssetFactory<Dynamics, Config&>::GetRegistry() {
    return Registry();
  }

  // Explicit template instantiation to ensure single registry across library boundaries
  template class AssetFactory<Dynamics, Config&>;

}  // namespace lupnt
