/**
 * @file state.h
 * @author Stanford NAV LAB
 * @brief  SPICE Interface functions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once
#include <spdlog/spdlog.h>

#include "lupnt/core/definitions.h"
#include "lupnt/dynamics/dynamics.h"
#include "lupnt/filters/filter.h"
#include "lupnt/states/state.h"

namespace lupnt {

  // define Params
  enum EstType { ESTIMATED, CONSIDERED, FIXED };

  /// @brief Concatenation of multiple `State`/`Dynamics`/`ParamState` triples
  /// (e.g. one per agent or per orbit+clock sub-state) into a single combined
  /// state vector, dynamics function, and process-noise function for use by a
  /// `Filter` (EKF/UKF/batch least squares).
  ///
  /// Built up incrementally via `Add()` (one call per sub-state), then queried
  /// via `GetStmState`/`GetDynamicsFunction`/`GetProcessNoiseFunction` to
  /// obtain the propagation function and initial state/STM the filter
  /// integrates each step. Per-parameter `EstType` (ESTIMATED / CONSIDERED /
  /// FIXED) and first-order Gauss-Markov time constants (`tauq`) control which
  /// parameters appear in the filter's state-transition matrix and how their
  /// process noise evolves.
  class JointState {
  private:
    std::vector<State> states_;
    std::vector<Ptr<Dynamics>> dynamics_;
    std::vector<ParamState> params_;
    std::vector<std::vector<EstType>> est_types_;
    std::vector<Ptr<ProcessNoiseFunction>> proc_noise_;
    MatXd est_matrix_;  // K x N matrix where M is the number of estimated parameters and N is the
                        // number of states (N = K + M)
    MatXd consider_matrix_;  // M x N matrix where M is the number of consider parameters and N is
                             // the number of states
    std::vector<std::vector<std::pair<double, double>>>
        tauqs_;  // Time constants and sigmas for first-order Gauss Markov

    double tau_default_ = 1e12;  // Default time constant for first-order Gauss-Markov

    int size_ = 0;
    int state_size_ = 0;
    int param_size_ = 0;
    int considered_param_size_ = 0;
    int estimated_param_size_ = 0;
    int fixed_param_size_ = 0;
    int stm_size_ = 0;

  public:
    int GetSize() const { return state_size_ + param_size_; };
    int GetStateSize() const { return state_size_; };
    int GetParamSize() const { return param_size_; };
    int GetStmSize() const { return stm_size_; };
    int GetEstimatedParamSize() const { return estimated_param_size_; };
    int GetConsideredParamSize() const { return considered_param_size_; };
    int GetFixedParamSize() const { return fixed_param_size_; };
    State GetState() const;
    State GetStmState() const;
    ParamState GetParams() const;

    /**
     * @brief Add a new state to the joint state
     * @param state The state to add
     * @param dynamics The dynamics associated with the state
     * @param proc_noise The process noise function associated with the state
     * @param param The parameters associated with the state
     * @param est_types The estimation types for the parameters
     * @param tauq The time constants and sigmas for the parameters (first-order Gauss-Markov)
     */
    void Add(const State& state, Ptr<Dynamics>&& dynamics,
             Ptr<ProcessNoiseFunction>&& proc_noise = nullptr,
             const ParamState& param = ParamState(0), std::vector<EstType> est_types = {},
             std::vector<std::pair<double, double>> tauq = {});

    // void Add(State state, Ptr<DynamicsWithParams> dynamics,
    //                          Ptr<FilterProcessNoiseFunction> proc_noise = nullptr);

    FilterDynamicsFunction GetDynamicsFunction();
    ProcessNoiseFunction GetProcessNoiseFunction();
    MatXd GetProcessNoiseMappingMatrix();
  };

}  // namespace lupnt
