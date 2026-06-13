#pragma once

#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>
#include <functional>
#include <vector>

#include "lupnt/core/config.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  /**
   * @brief Configuration for batch filter algorithm
   */
  struct BatchFilterConfig {
    bool use_weights = true;
    bool use_initialization = true;
    double convergence_tol = 1e-9;
    int max_iterations = 1000;

    // Constructor from YAML
    BatchFilterConfig() = default;
    BatchFilterConfig(const Config& config) {
      if (config["use_weights"]) use_weights = config["use_weights"].as<bool>();
      if (config["use_initialization"])
        use_initialization = config["use_initialization"].as<bool>();
      if (config["convergence_tol"]) convergence_tol = config["convergence_tol"].as<double>();
      if (config["max_iterations"]) max_iterations = config["max_iterations"].as<int>();
    }
  };

  /**
   * @brief Results from batch filter estimation
   */
  struct BatchFilterResults {
    VecXd state_estimate;    // Final state estimate
    MatXd state_covariance;  // Final state covariance
    VecXd state_errors;      // State error magnitudes (if true state provided)
    double pdop;             // Position Dilution of Precision
    int iterations;          // Number of iterations to convergence
    bool converged;          // Whether filter converged

    BatchFilterResults() : pdop(0.0), iterations(0), converged(false) {}
  };

  /**
   * @brief Measurement model function type
   *
   * Takes current state estimate and measurement index, returns predicted measurement and Jacobian
   *
   * @param state Current state estimate [n_state]
   * @param meas_idx Index of measurement to predict
   * @return std::pair<VecXd, MatXd> (predicted_measurement [n_meas], jacobian [n_meas x n_state])
   */
  using MeasurementModelFunction = std::function<std::pair<VecXd, MatXd>(const VecXd&, int)>;

  /**
   * @brief Generic batch filter using iterative weighted least squares
   *
   * Solves for state that best fits all measurements using iterative nonlinear least squares.
   * Works with any measurement model provided via function pointer.
   *
   * @param initial_state Initial state estimate [n_state]
   * @param initial_covariance Initial state covariance diagonal [n_state] (for initialization
   * constraints)
   * @param measurements Vector of measurements [n_measurements] each of size [n_meas_per_obs]
   * @param measurement_weights Vector of measurement weights [n_measurements] each of size
   * [n_meas_per_obs] (1/sigma^2)
   * @param measurement_model Function that predicts measurements and computes Jacobians
   * @param config Batch filter configuration
   * @param true_state Optional true state for error computation [n_state]
   * @return BatchFilterResults
   */
  BatchFilterResults RunBatchFilter(const VecXd& initial_state, const VecXd& initial_covariance,
                                    const std::vector<VecXd>& measurements,
                                    const std::vector<VecXd>& measurement_weights,
                                    MeasurementModelFunction measurement_model,
                                    const BatchFilterConfig& config,
                                    const VecXd& true_state = VecXd());

  // ============================================================================
  // Helper Functions
  // ============================================================================

  /**
   * @brief Solve weighted least squares problem
   *
   * @param H Design matrix [n_obs x n_state]
   * @param residuals Measurement residuals [n_obs]
   * @param weights Measurement weights [n_obs] (1/sigma^2)
   * @return VecXd State correction [n_state]
   */
  VecXd SolveWeightedLeastSquares(const MatXd& H, const VecXd& residuals, const VecXd& weights);

  /**
   * @brief Calculate Position Dilution of Precision (PDOP)
   *
   * Assumes first 3 state elements are position
   *
   * @param covariance_matrix State covariance matrix [n_state x n_state]
   * @return double PDOP value (sqrt of trace of position block)
   */
  double CalculatePDOP(const MatXd& covariance_matrix);

  /**
   * @brief Add initialization constraints to weighted least squares problem
   *
   * Adds soft constraints to keep state estimate near initial guess
   *
   * @param H_meas Measurement design matrix [n_meas x n_state]
   * @param residuals_meas Measurement residuals [n_meas]
   * @param weights_meas Measurement weights [n_meas]
   * @param initial_state Initial state estimate [n_state]
   * @param current_state Current state estimate [n_state]
   * @param init_weights Initialization weights [n_state] (1/sigma^2)
   * @return std::tuple<MatXd, VecXd, VecXd> (H_augmented, residuals_augmented, weights_augmented)
   */
  std::tuple<MatXd, VecXd, VecXd> AddInitializationConstraints(
      const MatXd& H_meas, const VecXd& residuals_meas, const VecXd& weights_meas,
      const VecXd& initial_state, const VecXd& current_state, const VecXd& init_weights);

}  // namespace lupnt
