/**
 * @file StateEstimationApp.h
 * @author Stanford NAV LAB
 * @brief Base class for State Estimation Application
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once
#include <lupnt/agents/agent.h>
#include <lupnt/agents/application.h>
#include <lupnt/agents/comm_device.h>
#include <lupnt/core/file.h>
#include <lupnt/measurements/gnss_measurement.h>
#include <lupnt/measurements/gnss_receiver.h>
#include <lupnt/numerics/filters.h>
#include <lupnt/physics/state.h>

#include <memory>

namespace lupnt {

class StateEstimationApp : public Application {
 private:
  double epoch0_;  // start epoch in TAU
  double epoch_;   // curent epoch
  double t_;       // Current time [s]

  std::shared_ptr<Agent> agent_;
  std::shared_ptr<IFilter> filter_;
  std::vector<IDynamics *> dynamics_vec_;
  DynamicsFunction dynamics_func_;
  MeasurementFunction meas_func_;
  JointState state_vec_;

  std::shared_ptr<DataHistory> data_history_;

 public:
  void SetAgent(std::shared_ptr<Agent> agent) { this->agent_ = agent; }
  void SetDataHistory(std::shared_ptr<DataHistory> data_history) {
    this->data_history_ = data_history;
  };

  double GetInitialEpoch() { return epoch0_; };
  double GetCurrentEpoch() { return epoch_; };
  double GetCurrrentTime() { return t_; };

  void SetDynamicsFunction();

  void Setup();
  void Step(double t_end);  // execute to step t
};
};  // namespace lupnt
