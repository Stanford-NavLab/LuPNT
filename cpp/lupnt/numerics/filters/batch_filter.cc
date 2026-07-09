#include "lupnt/numerics/filters/batch_filter.h"

#include <Eigen/QR>
#include <cmath>
#include <iostream>

namespace lupnt {

  BatchFilterResults RunBatchFilter(const VecXd& initial_state, const VecXd& initial_covariance,
                                    const std::vector<VecXd>& measurements,
                                    const std::vector<VecXd>& measurement_weights,
                                    MeasurementModelFunction measurement_model,
                                    const BatchFilterConfig& config, const VecXd& true_state) {
    BatchFilterResults results;

    // Initialize state estimate
    VecXd state_estimate = initial_state;
    const int n_state = initial_state.size();
    const int n_measurements = measurements.size();

    if (n_measurements == 0) {
      results.state_estimate = initial_state;
      results.state_covariance = initial_covariance.asDiagonal();
      results.converged = false;
      results.pdop = std::numeric_limits<double>::infinity();
      return results;
    }

    // Count total number of measurement elements
    int total_meas_size = 0;
    for (const auto& meas : measurements) {
      total_meas_size += meas.size();
    }

    // Iterative nonlinear least squares
    bool converged = false;
    int iteration = 0;

    for (iteration = 0; iteration < config.max_iterations; iteration++) {
      // Build measurement vectors and design matrix
      VecXd all_residuals(total_meas_size);
      VecXd all_weights(total_meas_size);
      MatXd all_jacobians(total_meas_size, n_state);

      int meas_idx = 0;
      int elem_idx = 0;

      for (const auto& measurement : measurements) {
        // Get predicted measurement and Jacobian from model
        auto [predicted_meas, jacobian] = measurement_model(state_estimate, meas_idx);

        // Check dimensions
        if (predicted_meas.size() != measurement.size()) {
          throw std::runtime_error("Measurement size mismatch in batch filter");
        }
        if (jacobian.rows() != measurement.size() || jacobian.cols() != n_state) {
          throw std::runtime_error("Jacobian size mismatch in batch filter");
        }

        // Compute residuals: observed - predicted
        VecXd residual = measurement - predicted_meas;

        // Store in combined vectors
        int meas_size = measurement.size();
        all_residuals.segment(elem_idx, meas_size) = residual;
        all_weights.segment(elem_idx, meas_size) = measurement_weights[meas_idx];
        all_jacobians.block(elem_idx, 0, meas_size, n_state) = jacobian;

        elem_idx += meas_size;
        meas_idx++;
      }

      // Record pre-correction diagnostics for this iteration (measurement residuals
      // only, matching the "raw" fit quality before any initialization pull is added).
      const double weighted_rms
          = std::sqrt((all_weights.array() * all_residuals.array().square()).mean());

      // Add initialization constraints if enabled
      MatXd H_final = all_jacobians;
      VecXd residuals_final = all_residuals;
      VecXd weights_final = all_weights;

      if (config.use_initialization) {
        VecXd init_weights = 1.0 / initial_covariance.array();  // 1/sigma^2
        auto [H_aug, res_aug, w_aug] = AddInitializationConstraints(
            all_jacobians, all_residuals, all_weights, initial_state, state_estimate, init_weights);
        H_final = H_aug;
        residuals_final = res_aug;
        weights_final = w_aug;
      }

      // Solve weighted least squares
      VecXd state_correction;
      if (config.use_weights) {
        state_correction = SolveWeightedLeastSquares(H_final, residuals_final, weights_final);
      } else {
        // Unweighted least squares
        state_correction
            = (H_final.transpose() * H_final).ldlt().solve(H_final.transpose() * residuals_final);
      }

      results.iteration_history.push_back(
          {iteration, state_estimate, state_correction.norm(), weighted_rms});

      // Update state estimate
      state_estimate += state_correction;

      // Check convergence
      if (state_correction.norm() < config.convergence_tol) {
        converged = true;
        break;
      }
    }

    // Compute final covariance (measurement Jacobians only, not initialization)
    MatXd final_jacobians(total_meas_size, n_state);
    VecXd final_weights(total_meas_size);

    int elem_idx = 0;
    for (int meas_idx = 0; meas_idx < n_measurements; meas_idx++) {
      auto [predicted_meas, jacobian] = measurement_model(state_estimate, meas_idx);
      int meas_size = measurements[meas_idx].size();

      final_jacobians.block(elem_idx, 0, meas_size, n_state) = jacobian;
      final_weights.segment(elem_idx, meas_size) = measurement_weights[meas_idx];

      elem_idx += meas_size;
    }

    // Compute covariance matrix
    MatXd information_matrix;
    if (config.use_weights) {
      MatXd W = final_weights.asDiagonal();
      information_matrix = final_jacobians.transpose() * W * final_jacobians;
    } else {
      information_matrix = final_jacobians.transpose() * final_jacobians;
    }

    // Compute covariance (pseudo-inverse for rank-deficient cases)
    MatXd covariance_matrix;
    if (information_matrix.rows() > 0) {
      Eigen::JacobiSVD<MatXd> svd(information_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
      covariance_matrix = svd.solve(MatXd::Identity(n_state, n_state));
    } else {
      covariance_matrix = MatXd::Identity(n_state, n_state) * 1e6;  // Large uncertainty
    }

    // Fill results
    results.state_estimate = state_estimate;
    results.state_covariance = covariance_matrix;
    results.iterations = iteration + 1;
    results.converged = converged;
    results.pdop = CalculatePDOP(covariance_matrix);

    // Compute state errors if true state provided
    if (true_state.size() > 0 && true_state.size() == n_state) {
      results.state_errors = (state_estimate - true_state).cwiseAbs();
    }

    return results;
  }

  VecXd SolveWeightedLeastSquares(const MatXd& H, const VecXd& residuals, const VecXd& weights) {
    // Weighted least squares: solve (H^T W H) x = H^T W residuals
    MatXd W = weights.asDiagonal();
    MatXd HTW = H.transpose() * W;
    MatXd HTWH = HTW * H;
    VecXd HTWr = HTW * residuals;

    // Use LDLT decomposition for numerical stability
    return HTWH.ldlt().solve(HTWr);
  }

  double CalculatePDOP(const MatXd& covariance_matrix) {
    if (covariance_matrix.rows() < 3 || covariance_matrix.cols() < 3) {
      return std::numeric_limits<double>::infinity();
    }

    // PDOP = sqrt(trace of position covariance block)
    // Assumes first 3 state elements are position
    Mat3d pos_cov = covariance_matrix.block<3, 3>(0, 0);
    return std::sqrt(pos_cov.trace());
  }

  std::tuple<MatXd, VecXd, VecXd> AddInitializationConstraints(
      const MatXd& H_meas, const VecXd& residuals_meas, const VecXd& weights_meas,
      const VecXd& initial_state, const VecXd& current_state, const VecXd& init_weights) {
    const int n_meas = H_meas.rows();
    const int n_state = H_meas.cols();

    // Create identity matrix for initialization constraints
    MatXd H_init = MatXd::Identity(n_state, n_state);
    VecXd residuals_init = initial_state - current_state;  // Pull toward initial state

    // Combine measurement and initialization constraints
    MatXd H_combined(n_meas + n_state, n_state);
    H_combined.topRows(n_meas) = H_meas;
    H_combined.bottomRows(n_state) = H_init;

    VecXd residuals_combined(n_meas + n_state);
    residuals_combined.head(n_meas) = residuals_meas;
    residuals_combined.tail(n_state) = residuals_init;

    VecXd weights_combined(n_meas + n_state);
    weights_combined.head(n_meas) = weights_meas;
    weights_combined.tail(n_state) = init_weights;

    return std::make_tuple(H_combined, residuals_combined, weights_combined);
  }

}  // namespace lupnt
