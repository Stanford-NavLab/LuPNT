/**
 * @file Filters.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2023-09-09
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

// Autodiff includes
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

// Eigen includes
#include <Eigen/Dense>
#include <Eigen/QR>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

namespace ad = autodiff;
namespace LPT {

// Dynamics and Measurement Function

/**
 * @brief Dynamics function
 *
 * @param x State
 * @param t_curr Current time
 * @param t_end End time
 * @param Phi STM of the dynamics
 */
typedef std::function<ad::VectorXreal(const ad::VectorXreal, ad::real t_curr,
                                      ad::real t_end, Eigen::MatrixXd &)>
    DynamicsFunction;

/**
 * @brief Process noise function
 * 
 * @param x State
 * @param t_curr Current time
 * @param t_end End time
 * @return Eigen::MatrixXd Process noise covariance
 * 
 */
typedef std::function<Eigen::MatrixXd(const ad::VectorXreal, ad::real t_curr,
                                      ad::real t_end)>
    ProcessNoiseFunction;

/**
 * @brief Measurement function
 * 
 * @param x State
 * @param H Measurement matrix
 * @param R Measurement noise covariance
 * 
 */
typedef std::function<ad::VectorXreal(const ad::VectorXreal, Eigen::MatrixXd &,
                                      Eigen::MatrixXd &)>
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
  ad::VectorXreal x;     // Updated state
  ad::real t_curr;       // Current time
  ad::VectorXreal xbar;  // Predicted state
  ad::VectorXreal dx;    // State update

  Eigen::MatrixXd P;     // Updated state cov
  Eigen::MatrixXd Pbar;  // Predicted state cov

  // Eigen::MatrixXd F;  // Jacobian of the dynamics
  Eigen::MatrixXd H;  // Measurement matrix
  Eigen::MatrixXd S;  // Innovation cov
  Eigen::MatrixXd K;  // Kalman gain
  Eigen::MatrixXd I;  // Identity matrix

  void Initialize(const ad::VectorXreal &x0, const Eigen::MatrixXd &P0) {
    x = x0;
    P = P0;
  }

  ad::VectorXreal GetPredictedStateEstimate(Eigen::MatrixXd &Pbar) {
    Pbar = this->Pbar;
    return this->xbar;
  }

  ad::VectorXreal GetUpdatedStateEstimate(Eigen::MatrixXd &P) {
    Pbar = this->P;
    return this->x;
  }

  void Predict(ad::real t_end);
  void Update(ad::VectorXreal z_obs);

  void Step(ad::real t_end, ad::VectorXreal z_obs) {
    Predict(t_end);
    Update(z_obs);
  }
};

}  // namespace LPT