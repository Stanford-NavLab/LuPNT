/**
 * @file GnssStateEstimationApp.h
 * @author Stanford NAV LAB
 * @brief State Estimation Application Using Gnss Measurement
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

#include <memory>

namespace lupnt {

class GnssStateEstimationApp : public Application {
 private:
  double epoch0;     // start epoch in TAU
  double epoch;      // curent e
  double t = 0.0;    // Current time [s]
  double dt = 1.0;   // Integration time step [s]
  double Dt = 10.0;  // Propagation time step [s]

  double frequency = 10.0;

  // Estimation
  int n_state = 8;
  double pos_err = 1.0;            // Position error [km]
  double vel_err = 1e-3;           // Velocity error [km/s]
  double clk_bias_err = 1e-6;      // Clock bias error [s]
  double clk_drift_err = 1e-9;     // Clock drift error [s/s]
  double sigma_acc = 2e-6;         // Process noise [km/s^2]
  double sigma_range = 5e-3;       // Range measurement noise [km]
  double sigma_range_rate = 1e-6;  // Range rate measurement noise [km/s]

  ad::Vector6real rv_est;
  ad::Vector2real clk_est;
  ad::Vector6real rv_pred_only;
  ad::Vector2real clk_pred_only;

  Eigen::MatrixXd P;
  Eigen::MatrixXd P_pred_only;

  Eigen::Matrix6d Q_rv;
  Eigen::Matrix2d Q_clk;
  Eigen::MatrixXd Q;

  // State transition matrices
  Eigen::Matrix6d Phi_rv;
  Eigen::Matrix6d Phi_rv_pred_only;
  Eigen::Matrix2d Phi_clk;
  Eigen::Matrix2d Phi_clk_pred_only;

  std::shared_ptr<Agent> agent;
  std::shared_ptr<IOrbitDynamics> dyn;
  std::shared_ptr<GnssReceiver> receiver;
  std::shared_ptr<IFilter> filter;
  MeasurementFunction meas_func;

  ad::VectorXreal x_est;
  ad::VectorXreal P_est;

  std::shared_ptr<DataHistory> data_history;

 public:
  void SetAgent(std::shared_ptr<Agent> agent) { this->agent = agent; }
  void SetDynamics(std::shared_ptr<IOrbitDynamics> dyn) { this->dyn = dyn; }
  void SetReceiver(std::shared_ptr<GnssReceiver> receiver) {
    this->receiver = receiver;
  }
  void SetDataHistory(std::shared_ptr<DataHistory> data_history) {
    this->data_history = data_history;
  }

  void SetFrequency(double frequency) { this->frequency = frequency; }
  double GetFrequency() { return frequency; }

  void Setup();
  void Step(double t);
};
};  // namespace lupnt
