#include "lupnt/numerics/filters/filter.h"

#include <iostream>

#include "lupnt/core/data_logger.h"

namespace lupnt {

  void Filter::Log(Real time) {
    (void)time;
    DataLogger::Log(fmt::format("{}/time", name_), t_);
    // for (int i = 0; i < x_.size(); i++) {
    //   auto name = x_.GetNames()[i];
    //   auto unit = x_.GetUnits()[i];
    //   std::replace(unit.begin(), unit.end(), '/', '_');
    //   DataLogger::Log(fmt::format("{}/state/{}_{}", name_, name, unit), x_(i));
    //   DataLogger::Log(fmt::format("{}/state_prior/{}_{}", name_, name, unit), x_prior_(i));
    //   DataLogger::Log(fmt::format("{}/state_post/{}_{}", name_, name, unit), x_post_(i));
    //   DataLogger::Log(fmt::format("{}/sigmas/{}_{}", name_, name, unit), P_(i, i));
    //   DataLogger::Log(fmt::format("{}/sigmas_prior/{}_{}", name_, name, unit), P_prior_(i, i));
    //   DataLogger::Log(fmt::format("{}/sigmas_post/{}_{}", name_, name, unit), P_post_(i, i));
    // }
    DataLogger::Log(fmt::format("{}/state", name_), x_);
    DataLogger::Log(fmt::format("{}/state_prior", name_), x_prior_);
    DataLogger::Log(fmt::format("{}/state_post", name_), x_post_);
    DataLogger::Log(fmt::format("{}/sigmas", name_), P_.diagonal());
    DataLogger::Log(fmt::format("{}/sigmas_prior", name_), P_prior_.diagonal());
    DataLogger::Log(fmt::format("{}/sigmas_post", name_), P_post_.diagonal());
  };

  double Filter::ComputeResidualRMS(VecXd z_true, const State& x_est) {
    MatXd H;
    MatXd R;
    if (z_true.size() == 0) {
      return 0.0;
    }
    VecXd z_pred = f_meas_(x_est, &H, &R);
    VecXd dz = z_pred - z_true;
    VecXd R_diag = R.diagonal();
    double sum_dz = 0.0;
    for (int i = 0; i < dz.size(); i++) {
      // if dz is NaN, skip
      if (std::isnan(dz(i)) || R_diag(i) <= 0.0) {
        continue;
      }
      sum_dz += (dz(i) * dz(i)) / R_diag(i);
    }
    // std::cout << "dz: " << dz.transpose() << std::endl;
    // std::cout << "R diag: " << R_diag.transpose() << std::endl;
    // std::cout << "Residual RMS: " << sum_dz << std::endl;
    return sum_dz;
  };

}  // namespace lupnt
