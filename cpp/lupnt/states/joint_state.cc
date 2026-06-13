
#include "lupnt/states/joint_state.h"
namespace lupnt {

  State JointState::GetState() const {
    State x(state_size_);
    int idx = 0;
    for (auto& state : states_) {
      for (int i = 0; i < state.size(); i++) {
        x(idx) = state(i);
        x.GetUnits()[idx] = state.GetUnits()[i];
        x.GetNames()[idx] = state.GetNames()[i];
        idx++;
      }
    }
    return x;
  }

  /**
   * \returns the concatenated parameter state
   */
  ParamState JointState::GetParams() const {
    ParamState p(param_size_);
    int idx = 0;
    for (auto& param : params_) {
      if (param.size() == 0) {
        continue;
      }
      for (int i = 0; i < param.size(); i++) {
        p(idx) = param(i);
        p.GetUnits()[idx] = param.GetUnits()[i];
        p.GetNames()[idx] = param.GetNames()[i];
        idx++;
      }
    }
    return p;
  }

  State JointState::GetStmState() const {
    State x(stm_size_);
    int idx = 0;
    // Add states
    for (auto& state : states_) {
      for (int i = 0; i < state.size(); i++) {
        x(idx) = state(i);
        x.GetUnits()[idx] = state.GetUnits()[i];
        x.GetNames()[idx] = state.GetNames()[i];
        idx++;
      }
    }
    // Add estimated and considered parameters
    for (size_t j = 0; j < params_.size(); j++) {
      const auto& param = params_[j];
      const auto& est_types = est_types_[j];
      for (size_t k = 0; k < param.size(); k++) {
        if (est_types[k] == ESTIMATED || est_types[k] == CONSIDERED) {
          x(idx) = param(k);
          x.GetUnits()[idx] = param.GetUnits()[k];
          x.GetNames()[idx] = param.GetNames()[k];
          idx++;
        }
      }
    }
    return x;
  }

  void JointState::Add(const State& state, Ptr<Dynamics>&& dynamics,
                       Ptr<ProcessNoiseFunction>&& proc_noise, const ParamState& param,
                       std::vector<EstType> est_types,
                       std::vector<std::pair<double, double>> tauq) {
    states_.push_back(state);
    dynamics_.emplace_back(std::move(dynamics));
    if (proc_noise)
      proc_noise_.push_back(std::move(proc_noise));
    else
      proc_noise_.push_back(nullptr);
    params_.emplace_back(param);

    // Parameters
    if (param.size() > 1 && est_types.size() == 1) {
      // If only one est_type is provided, apply it to all parameters
      spdlog::info(
          "Single estimation type provided for multiple parameters. Applying to all parameters.");
      est_types = std::vector<EstType>(param.size(), est_types[0]);
    } else if (est_types.size() != param.size()) {
      spdlog::warn(
          "Estimation type size does not match parameter size. All parameters set to FIXED.");
      est_types = std::vector<EstType>(param.size(), FIXED);
    }
    est_types_.emplace_back(est_types);

    // First-Order Gauss-Markov Process Noise
    if (tauq.size() == 1 && param.size() > 1) {
      // If only one tau is provided, apply it to all parameters
      spdlog::info(
          "Single time constant provided for multiple parameters. Applying to all parameters.");
      tauq = std::vector<std::pair<double, double>>(param.size(),
                                                    std::make_pair(tauq[0].first, tauq[0].second));
    } else if (tauq.size() != param.size()) {
      spdlog::warn(
          "Time constant size does not match parameter size. All time constants set to 1e10 "
          "(infinity).");
      tauq = std::vector<std::pair<double, double>>(
          param.size(), std::make_pair(tau_default_, 0.0));  // Set to large value
    }
    tauqs_.emplace_back(tauq);
    state_size_ += state.size();
    param_size_ += param.size();
    considered_param_size_ += std::count(est_types.begin(), est_types.end(), CONSIDERED);
    estimated_param_size_ += std::count(est_types.begin(), est_types.end(), ESTIMATED);
    fixed_param_size_ += std::count(est_types.begin(), est_types.end(), FIXED);
    stm_size_ = state_size_ + estimated_param_size_ + considered_param_size_;
  }

  FilterDynamicsFunction JointState::GetDynamicsFunction() {
    // State propagation function
    // State = [state1, state2, ..., param1, param2, ...]
    // where param are only estimated and considered parameters
    // STM = d(State + param)/d(State + param)
    // The other parameters are treated as constants (taken from params_)
    FilterDynamicsFunction func = [this](const State& x0, Real t0, Real tf, const State* u,
                                         MatXd* Phi) {
      int i = 0;    // index for states_
      int idx = 0;  // index for state vector
      int param_idx = 0;
      State xf = x0;
      Real dt = tf - t0;

      // Resize STM
      if (Phi != nullptr) {
        Phi->resize(stm_size_, stm_size_);
        // Set Eye matrix
        Phi->setIdentity();
      }

      // States propagation
      for (auto& state : states_) {
        int n = state.size();
        State x_tmp = x0.segment(idx, n);

        // Insert Parameters to State
        ParamState param = params_[i];
        int m = 0;  // number of estimated and considered parameters for this dynamics
        int jj = 0;
        for (int k = 0; k < param.size(); k++) {
          // Set only estimated and considered parameters
          if (est_types_[i][k] == ESTIMATED || est_types_[i][k] == CONSIDERED) {
            // Constant Dynamics
            double tau = tauqs_[i][k].first;
            if (tau < tau_default_) {
              // First-order Gauss-Markov process
              param(k) = exp(-dt / tau) * x0(state_size_ + param_idx + jj);
            } else {
              // Constant parameter
              param(k) = x0(state_size_ + param_idx + jj);
            }
            params_[i](k) = param(
                k);  // Update the parameter in params_ with the current value from state vector
            m++;
            jj++;
          }
        }

        if (Phi != nullptr) {
          MatXd Phi_state(n, n);
          MatXd Phi_param(n, param_size_);  // over-allocated, will use only m columns
          if (m > 0) {
            // State and Parameter propagation
            x_tmp = dynamics_[i]->PropagateWithParams(x_tmp, t0, tf, param, u, &Phi_state,
                                                      &Phi_param);

            // Param->State Mapping (Only estimated and considered parameters) and Param->Param
            // Mapping
            jj = 0;
            for (int k = 0; k < param.size(); k++) {  // loop over all parameters (including fixed)
              if (est_types_[i][k] == FIXED) {
                continue;
              }
              // Param->State Mapping
              Phi->block(0, state_size_ + param_idx + jj, n, 1) = Phi_param.col(k);

              // Param->Param Mapping
              if (tauqs_[i][k].first < tau_default_) {
                // First-order Gauss-Markov process
                (*Phi)(state_size_ + param_idx + jj, state_size_ + param_idx + jj)
                    = exp(-dt / tauqs_[i][k].first);
              } else {
                // Constant parameter
                (*Phi)(state_size_ + param_idx + jj, state_size_ + param_idx + jj) = 1.0;
              }

              // Increment
              jj++;
            }

          } else {
            x_tmp = dynamics_[i]->Propagate(x_tmp, t0, tf, u, &Phi_state);
          }
          // State->State Mapping
          Phi->block(idx, idx, n, n) = Phi_state;

        } else {
          x_tmp = dynamics_[i]->Propagate(x_tmp, t0, tf, u);
        }
        xf.segment(idx, n) = x_tmp;
        idx += n;
        param_idx += m;
        i++;
      }
      return xf;
    };
    return func;
  };

  ProcessNoiseFunction JointState::GetProcessNoiseFunction() {
    ProcessNoiseFunction func = [this](const State& x, Real t0, Real tf) {
      MatXd Q = MatXd::Zero(stm_size_, stm_size_);
      int idx = 0;
      int i = 0;
      int param_idx = 0;

      for (auto& state : states_) {
        // States
        int n = state.size();
        if (proc_noise_[idx]) {
          State x_tmp = x.segment(idx, n);
          MatXd Q_tmp = (*proc_noise_[i])(x_tmp, t0, tf);
          Q.block(idx, idx, n, n) = Q_tmp;
        }

        // Parameters
        int m = params_[i].size();
        int jj = 0;
        if (m > 0) {
          MatXd Q_param = MatXd::Zero(m, m);
          for (int k = 0; k < m; k++) {
            if (est_types_[i][k] == ESTIMATED || est_types_[i][k] == CONSIDERED) {
              double tau = tauqs_[i][k].first;
              double q = tauqs_[i][k].second;
              if (tau < tau_default_) {
                // First-order Gauss-Markov process
                double exp_factor = exp(-2 * (tf - t0) / tau);
                double Q_next = q * tau / 2 * (1 - exp_factor);
                Q(state_size_ + param_idx + jj, state_size_ + param_idx + jj) = Q_next;
              } else {
                // Constant parameter
                Q(state_size_ + param_idx + jj, state_size_ + param_idx + jj) = 0.0;
              }
              jj++;
            }
          }
        }

        // Increment
        idx += n;
        i++;
        param_idx += m;
      }
      return Q;
    };
    return func;
  }

}  // namespace lupnt
