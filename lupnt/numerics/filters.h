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

// Autodiff includes



// Eigen includes
#include <lupnt/core/constants.h>


#include <Eigen/QR>




namespace lupnt {

// Dynamics and Measurement Function

/**
 * @brief Dynamics function
 *
 * @param x State
 * @param t_curr Current time
 * @param t_end End time
 * @param Phi STM of the dynamics
 */
typedef std::function<VectorX(const VectorX, real t_curr, real t_end,
                                  MatrixXd &)>
    DynamicsFunction;

/**
 * @brief Process noise function
 *
 * @param x State
 * @param t_curr Current time
 * @param t_end End time
 * @return MatrixXd Process noise covariance
 *
 */
typedef std::function<MatrixXd(const VectorX, real t_curr, real t_end)>
    ProcessNoiseFunction;

/**
 * @brief Measurement function
 *
 * @param x State
 * @param H Measurement matrix
 * @param R Measurement noise covariance
 *
 */
typedef std::function<VectorX(const VectorX, MatrixXd &, MatrixXd &)>
    MeasurementFunction;

class IFilter {
 public:
  virtual ~IFilter() = default;
  DynamicsFunction dynamics;
  ProcessNoiseFunction process_noise;
  MeasurementFunction measurement;

  void SetDynamicsFunction(DynamicsFunction dynamics) {
    this->dynamics = dynamics;
  }
  void SetProcessNoiseFunction(ProcessNoiseFunction process_noise) {
    this->process_noise = process_noise;
  }
  void SetMeasurementFunction(MeasurementFunction measurement) {
    this->measurement = measurement;
  }
};

/**
 * @brief Extended Kalman Filter
 *
 */
class EKF : public IFilter {
 public:
  VectorX x;     // Updated state
  real t_curr;       // Current time
  VectorX xbar;  // Predicted state
  VectorX dx;    // State update

  MatrixXd P;     // Updated state cov
  MatrixXd Pbar;  // Predicted state cov

  // MatrixXd F;  // Jacobian of the dynamics
  MatrixXd H;  // Measurement matrix
  MatrixXd S;  // Innovation cov
  MatrixXd K;  // Kalman gain
  MatrixXd I;  // Identity matrix

  void Initialize(const VectorX &x0, const MatrixXd &P0) {
    x = x0;
    P = P0;
  }

  VectorX GetPredictedStateEstimate(MatrixXd &Pbar) {
    Pbar = this->Pbar;
    return this->xbar;
  }

  VectorX GetUpdatedStateEstimate(MatrixXd &P) {
    Pbar = this->P;
    return this->x;
  }

  void Predict(real t_end);
  void Update(VectorX z_obs);

  void Step(real t_end, VectorX z_obs) {
    Predict(t_end);
    Update(z_obs);
  }
};

}  // namespace lupnt