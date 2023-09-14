/**
 * @file StateEstimationApp.h
 * @author
 * @brief
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once
#include <lupnt/agents/Agent.h>
#include <lupnt/agents/Application.h>
#include <lupnt/agents/CommDevice.h>
#include <lupnt/core/File.h>
#include <lupnt/measurements/GNSSMeasurement.h>
#include <lupnt/measurements/GNSSReceiver.h>
#include <lupnt/numerics/Filters.h>
#include <lupnt/physics/State.h>

#include <memory>

namespace LPT {

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
};  // namespace LPT
