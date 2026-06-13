#pragma once

#include "lupnt/core/definitions.h"

namespace lupnt {

  void PrintEKFProgressHeaderPVC();

  void PrintEKFProgressHeaderPVC(int n_sat);

  void PrintEKFProgressPVC(double t, double x_pos_err, double x_vel_err, double x_clk_bias_err);

  void PrintEKFProgressPVC(double t, const VecXd& est_err, int n_sat);

  void PrintEstimationStatistics(const VecXd& num_meas, const MatXd& error_mat, double ratio,
                                 int n_sat);
}  // namespace lupnt
