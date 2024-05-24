/**
 * @file filters.h
 * @author Stanford NAV LAB
 * @brief List of Filters
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <lupnt/core/constants.h>

#include <Eigen/QR>

namespace lupnt {

// Dynamics and Measurement Function

/**
 * @brief Dynamics function for the filter
 *
 * @param x State
 * @param t_curr Current time
 * @param t_end End time
 * @param Phi STM of the dynamics
 */
typedef std::function<VectorX(const VectorX, real t_curr, real t_end,
                              MatrixXd &)>
    FilterDynamicsFunction;

/**
 * @brief Process noise function for the filter
 *
 * @param x State
 * @param t_curr Current time
 * @param t_end End time
 * @return MatrixXd Process noise covariance
 *
 */
typedef std::function<MatrixXd(const VectorX, real t_curr, real t_end)>
    FilterProcessNoiseFunction;

/**
 * @brief Measurement function for the filter
 *
 * @param x State
 * @param H Measurement matrix
 * @param R Measurement noise covariance
 *
 */
typedef std::function<VectorX(const VectorX, MatrixXd &, MatrixXd &)>
    FilterMeasurementFunction;

class IFilter {
 public:
  virtual ~IFilter() = default;
  FilterDynamicsFunction dynamics_;
  FilterProcessNoiseFunction process_noise_;
  FilterMeasurementFunction measurement_;

  void SetDynamicsFunction(FilterDynamicsFunction dynamics) {
    dynamics_ = dynamics;
  }
  void SetProcessNoiseFunction(FilterProcessNoiseFunction process_noise) {
    process_noise_ = process_noise;
  }
  void SetMeasurementFunction(FilterMeasurementFunction measurement) {
    measurement_ = measurement;
  }
};

/**
 * @brief Extended Kalman Filter
 *
 */
class EKF : public IFilter {
 public:
  real t_curr_;     // Current time
  VectorX x_;       // Updated state
  VectorX xbar_;    // Predicted state
  VectorX dy_;      // Measurement residual
  VectorX dx_;      // State update
  VectorX z_true_;  // Observed measurement
  VectorX z_pred_;  // Predicted measurement

  MatrixXd P_;     // Updated state cov
  MatrixXd Pbar_;  // Predicted state cov
  MatrixXd Q_;     // Process noise cov

  MatrixXd H_;  // Measurement matrix
  MatrixXd S_;  // Innovation cov
  MatrixXd K_;  // Kalman gain
  MatrixXd R_;  // Measurement noise cov

  double outlier_threshold_ = 3.0;

  EKF() {}

  EKF(FilterDynamicsFunction dynamics, FilterProcessNoiseFunction process_noise,
      FilterMeasurementFunction measurement) {
    dynamics_ = dynamics;
    process_noise_ = process_noise;
    measurement_ = measurement;
  }

  void Initialize(const VectorX &x0, const MatrixXd &P0) {
    x_ = x0;
    P_ = P0;
  }

  VectorX GetPredictedStateEstimate(MatrixXd &Pbar) {
    Pbar = Pbar_;
    return xbar_;
  }
  VectorX GetPredictedStateEstimate() { return xbar_; }
  VectorX GetUpdatedStateEstimate(MatrixXd &Phat) {
    Phat = P_;
    return x_;
  }
  VectorX GetUpdatedStateEstimate() { return x_; }
  VectorX GetMeasurementResidual() { return dy_; }
  MatrixX GetKalmanGain() { return K_; }
  MatrixX GetMeasurementNoiseCov() { return R_; }
  MatrixX GetMeasurementJacobian() { return H_; }
  int GetMeasurementSize() { return H_.rows(); }

  void SetOutlierThreshold(double outlier_threshold) {
    assert(outlier_threshold >= 0 && "Outlier threshold must be positive");
    outlier_threshold_ = outlier_threshold;
  }

  /**
   * @brief Remove outliers from the measurement
   *
   * @param m   number of measurements
   * @param debug   debug flag
   * @return int   number of measurements after removing outliers
   */
  int RemoveOutliers(int m, bool debug = false);

  /**
   * @brief Predict the state to a given time
   *
   * @param t_end   end time
   */
  void Predict(real t_end);

  /**
   * @brief Update the state with a measurement
   *
   * @param z_obs   measurement
   * @param debug   debug flag
   */
  void Update(VectorX z_obs, bool debug = false);

  /**
   * @brief Update the state with a measurement
   *
   * @param t_end   end time
   * @param z_obs   measurement obtained at end time
   * @param debug   debug flag
   */
  void Step(real t_end, VectorX z_obs, bool debug = false);
};

}  // namespace lupnt