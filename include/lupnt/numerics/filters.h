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

#include <Eigen/QR>

#include "lupnt/core/constants.h"

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
  typedef std::function<VecX(const VecX, Real t_curr, Real t_end, MatXd &)> FilterDynamicsFunction;

  /**
   * @brief Process noise function for the filter
   *
   * @param x State
   * @param t_curr Current time
   * @param t_end End time
   * @return VecXd Process noise covariance
   *
   */
  typedef std::function<MatXd(const VecX, Real t_curr, Real t_end)> FilterProcessNoiseFunction;

  /**
   * @brief Measurement function for the filter
   *
   * @param x State
   * @param H Measurement matrix
   * @param R Measurement noise covariance
   *
   */
  typedef std::function<VecX(const VecX, MatXd &, MatXd &)> FilterMeasurementFunction;

  class IFilter {
  public:
    virtual ~IFilter() = default;
    FilterDynamicsFunction dynamics_;
    FilterProcessNoiseFunction process_noise_;
    FilterMeasurementFunction measurement_;

    void SetDynamicsFunction(FilterDynamicsFunction dynamics) { dynamics_ = dynamics; }
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
    Real t_curr_;  // Current time
    VecX x_;       // Updated state
    VecX xbar_;    // Predicted state
    MatXd Phi_;    // State transition matrix
    VecX dy_;      // Measurement residual
    VecX dx_;      // State update
    VecX z_true_;  // Observed measurement
    VecX z_pred_;  // Predicted measurement

    MatXd P_;     // Updated state cov
    MatXd Pbar_;  // Predicted state cov
    MatXd Q_;     // Process noise cov

    MatXd H_;  // Measurement matrix
    MatXd S_;  // Innovation cov
    MatXd K_;  // Kalman gain
    MatXd R_;  // Measurement noise cov

    double outlier_threshold_ = 3.0;

    EKF() {}

    EKF(FilterDynamicsFunction dynamics, FilterProcessNoiseFunction process_noise,
        FilterMeasurementFunction measurement) {
      dynamics_ = dynamics;
      process_noise_ = process_noise;
      measurement_ = measurement;
    }

    void Initialize(const VecX &x0, const MatXd &P0) {
      x_ = x0;
      P_ = P0;
    }

    VecX GetPredictedStateEstimate(MatXd &Pbar) {
      Pbar = Pbar_;
      return xbar_;
    }
    VecX GetPredictedStateEstimate() { return xbar_; }
    VecX GetUpdatedStateEstimate(MatXd &Phat) {
      Phat = P_;
      return x_;
    }
    VecX GetUpdatedStateEstimate() { return x_; }
    VecX GetMeasurementResidual() { return dy_; }
    MatX GetKalmanGain() { return K_; }
    MatX GetMeasurementNoiseCov() { return R_; }
    MatX GetMeasurementJacobian() { return H_; }
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
    void Predict(Real t_end);

    /**
     * @brief Update the state with a measurement
     *
     * @param z_obs   measurement
     * @param debug   debug flag
     */
    void Update(VecX z_obs, bool debug = false);

    /**
     * @brief Update the state with a measurement
     *
     * @param t_end   end time
     * @param z_obs   measurement obtained at end time
     * @param debug   debug flag
     */
    void Step(Real t_end, VecX z_obs, bool debug = false);
  };

}  // namespace lupnt
