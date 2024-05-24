/**
 * @file utils.cc
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2023-12-23
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

// lupnt includes
#include <lupnt/agents/agent.h>
#include <lupnt/agents/gnss_constellation.h>
#include <lupnt/agents/state_estimation_app.h>
#include <lupnt/core/file.h>
#include <lupnt/dynamics/dynamics.h>
#include <lupnt/measurements/gnss_channel.h>
#include <lupnt/measurements/gnss_measurement.h>
#include <lupnt/measurements/gnss_receiver.h>
#include <lupnt/measurements/transmission.h>
#include <lupnt/numerics/filters.h>
#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/clock.h>
#include <lupnt/physics/coord_converter.h>
#include <lupnt/physics/orbit_state.h>
#include <lupnt/physics/spice_interface.h>

// Eigen includes
#include <Eigen/QR>

// C++ includes
#include <matplot/matplot.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

/*******************************************
 * Functions for Filtering Examples
 *******************************************/
namespace lupnt {

static const int state_size = 8;

MatrixXd ConstructInitCovariance(double pos_err, double vel_err,
                                 double clk_bias_err, double clk_drift_err) {
  Matrix6d P_rv = Matrix6d::Zero();
  P_rv.block(0, 0, 3, 3) = Matrix3d::Identity() * pow(pos_err, 2);
  P_rv.block(3, 3, 3, 3) = Matrix3d::Identity() * pow(vel_err, 2);

  Matrix2d P_clk = Matrix2d::Zero();
  P_clk(0, 0) = pow(clk_bias_err, 2);
  P_clk(1, 1) = pow(clk_drift_err, 2);

  MatrixXd P0 = blkdiag(P_rv, P_clk);

  return P0;
};

void AddStateEstimationData(const std::shared_ptr<DataHistory> data_history,
                            const std::shared_ptr<Spacecraft> sat, EKF* ekf,
                            GnssConstellation* gps_const, GnssMeasurement* meas,
                            double t, double epoch) {
  // Navigation
  data_history->AddData("z_true", t, ekf->z_true_);
  data_history->AddData("z_pred", t, ekf->z_pred_);
  data_history->AddData("CN0", t, meas->GetCN0());

  data_history->AddData("vis_earth", t, meas->GetEarthOccultation());
  data_history->AddData("vis_moon", t, meas->GetMoonOccultation());
  data_history->AddData("vis_antenna", t, meas->GetMoonOccultation());
  data_history->AddData("vis_ionos", t, meas->GetMoonOccultation());

  // Moon spacecraft
  auto state = sat->GetCartesianGCRFStateAtEpoch(epoch);
  auto sate_mi = ConvertOrbitStateCoordSystem(state, epoch, Frame::MI);
  auto state_gcrf = ConvertOrbitStateCoordSystem(state, epoch, Frame::GCRF);
  data_history->AddData("rv_moon_mi", t, sate_mi->GetVector());
  data_history->AddData("rv_moon_gcrf", t, state_gcrf->GetVector());

  // Estimation
  data_history->AddData("rv", t, sat->GetOrbitState()->GetVector());
  data_history->AddData("rv_pred", t, ekf->xbar_.head(6));
  data_history->AddData("rv_est", t, ekf->x_.head(6));

  data_history->AddData("clk", t, sat->GetClockState().GetVector());
  data_history->AddData("clk_pred", t, ekf->xbar_.tail(2));
  data_history->AddData("clk_est", t, ekf->x_.tail(2));

  data_history->AddData("P_rv", t, ekf->P_.diagonal().segment(0, 6));
  data_history->AddData("P_clk", t, ekf->P_.diagonal().segment(6, 2));

  // GPS constellation
  for (int i = 0; i < gps_const->GetNumSatellites(); i++) {
    auto sate = gps_const->GetSatellite(i)->GetCartesianGCRFStateAtEpoch(epoch);
    auto sate_mi = ConvertOrbitStateCoordSystem(sate, epoch, Frame::MI);
    auto state_gcrf = ConvertOrbitStateCoordSystem(sate, epoch, Frame::GCRF);

    std::string name = "sat" + std::to_string(i);
    data_history->AddData(name + "_mi", t, sate->GetVector());
    data_history->AddData(name + "_gcrf", t, state_gcrf->GetVector());
  }

  // Bodies
  Vector6 vz6;
  vz6.setZero();
  data_history->AddData(
      "earth_mi", t,
      CoordConverter::Convert(epoch, vz6, Frame::GCRF, Frame::MI));
  data_history->AddData(
      "moon_gcrf", t,
      CoordConverter::Convert(epoch, vz6, Frame::MI, Frame::GCRF));
};

void PrintProgressHeader() {
  std::cout << "Run Simulation" << std::endl;
  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Time [min]  | Pos Err [m] | Vel Err [mm/s] | Clk Bias Err [ms]"
            << std::endl;
  std::cout << "--------------------------------------------------------------"
            << std::endl;
}

VectorXd ComputeEstimationErrors(const std::shared_ptr<Spacecraft> sat,
                                 EKF* ekf) {
  auto x_est = ekf->x_;
  auto x_true = sat->GetStateVector();

  double x_pos_err =
      1000 * (x_true.segment(0, 3) - x_est.segment(0, 3)).norm().val();
  double x_vel_err =
      1e6 * (x_true.segment(3, 3) - x_est.segment(3, 3)).norm().val();
  double x_clk_bias_err = 1e9 * abs((x_true(6) - x_est(6)).val());
  double x_clk_drift_err = 1e9 * abs((x_true(7) - x_est(7)).val());

  VectorXd est_err(4);
  est_err << x_pos_err, x_vel_err, x_clk_bias_err, x_clk_drift_err;

  return est_err;
}

void PrintProgress(double t, double x_pos_err, double x_vel_err,
                   double x_clk_bias_err) {
  std::cout.precision(5);
  std::cout << std::left << std::setw(12) << t / 60 << " " << std::left
            << std::setw(12) << x_pos_err << "  " << std::left << std::setw(14)
            << x_vel_err << "   " << std::left << std::setw(16)
            << x_clk_bias_err << std::endl;
};

void PrintEKFDebugInfo(EKF* ekf) {
  std::cout << "  Pbar : " << ekf->Pbar_.diagonal().transpose() << std::endl;
  std::cout << "  Q:  " << std::endl << ekf->Q_ << std::endl;
  std::cout << "  Kalman Gain: " << std::endl << ekf->K_ << std::endl;
  std::cout << "  H:  " << std::endl << ekf->H_ << std::endl;
  std::cout << "  R:  " << std::endl << ekf->R_ << std::endl;
  std::cout << "  S:  " << std::endl << ekf->S_ << std::endl;
  std::cout << "  Meas   Residuals: " << ekf->dy_.transpose() << std::endl;
  std::cout << "  Linear Residuals: "
            << (ekf->dy_ - ekf->H_ * ekf->dx_).transpose() << std::endl;
  std::cout << "  dx: " << ekf->dx_.transpose() << std::endl;
  std::cout << "  Phat: " << ekf->P_.diagonal().transpose() << std::endl;
  std::cout << "  " << std::endl;
}

/**
 * @brief Print Estimation Errors
 *
 * @param num_meas (n_time,)   Number of GPS measurements
 * @param error_mat (4, n_time)  Error Matrix (Position, Velocity, Clock Bias,
 * Clock Drift)
 */
void PrintEstimationStatistics(VectorXd num_meas, MatrixXd error_mat,
                               double data_ratio = 1.0) {
  int n_time = num_meas.size();
  Vector4d rms, means, stds, p68, p95, p99;

  if (error_mat.rows() != 4) {
    std::cout << "Wrong Matrix Size, Error matrix size must be (4 x timestep)"
              << std::endl;
    return;
  }

  // extract statistics range data
  int start_idx = (int)((1.0 - data_ratio) * n_time);
  int end_idx = n_time - 1;
  int n_range = end_idx - start_idx;

  VectorXd num_meas_range(n_range);
  MatrixXd error_mat_range(4, n_range);

  num_meas_range = num_meas.segment(start_idx, n_range);
  error_mat_range = error_mat.block(0, start_idx, 4, n_range);

  // compute statistics ----------------------------------------------------
  // rms
  for (int i = 0; i < 4; i++) {
    rms(i) = Rms(error_mat_range.row(i));
    means(i) = error_mat_range.row(i).mean();
    stds(i) = Std(error_mat_range.row(i));
    p68(i) = Percentile(error_mat_range.row(i), 0.68);
    p95(i) = Percentile(error_mat_range.row(i), 0.95);
    p99(i) = Percentile(error_mat_range.row(i), 0.99);
  }

  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "< Simulation Statistics (Last " << data_ratio * 100 << "%)>"
            << std::endl;
  std::cout << " " << std::endl;
  std::cout
      << "Statistics  | Position [m]  | Velocity [mm/s] | Clock Bias [ns] "
         "| Clk Drift [ns/s] "
      << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------------"
            << std::endl;

  std::cout.precision(5);
  std::cout << "RMS         | " << std::left << std::setw(16) << rms(0) << "  "
            << std::left << std::setw(16) << rms(1) << "   " << std::left
            << std::setw(16) << rms(2) << std::left << std::setw(16) << rms(3)
            << std::endl;
  std::cout << "Mean+-Std   | " << std::left << means(0) << "+-" << std::left
            << stds(0) << "  " << std::left << means(1) << "+-" << std::left
            << stds(1) << "   " << std::left << means(2) << "+-" << std::left
            << stds(2) << "   " << std::left << means(3) << "+-" << std::left
            << stds(3) << std::endl;
  std::cout << "68%         | " << std::left << std::setw(16) << p68(0) << "  "
            << std::left << std::setw(16) << p68(1) << "   " << std::left
            << std::setw(16) << p68(2) << std::left << std::setw(16) << p68(3)
            << std::endl;
  std::cout << "95%         | " << std::left << std::setw(16) << p95(0) << "  "
            << std::left << std::setw(16) << p95(1) << "   " << std::left
            << std::setw(16) << p95(2) << std::left << std::setw(16) << p95(3)
            << std::endl;
  std::cout << "99%         | " << std::left << std::setw(16) << p99(0) << "  "
            << std::left << std::setw(16) << p99(1) << "   " << std::left
            << std::setw(16) << p99(2) << std::left << std::setw(16) << p99(3)
            << std::endl;
  std::cout << "  " << std::endl;
}

void Plot3DTrajectory(const std::shared_ptr<DataHistory> data_history,
                      std::string state_type = "true") {
  using namespace matplot;

  // Plot trajectory
  auto fig = figure(true);
  auto ax = fig->current_axes();
  ax->hold(on);
  ax->grid(on);
  ax->xlabel("x [km]");
  ax->ylabel("y [km]");
  ax->zlabel("z [km]");

  std::vector<Timestamped<VectorXd>> data;
  if (state_type == "true") {
    ax->title("True Trajectory");
    data = data_history->GetData("rv");
  } else if (state_type == "est") {
    ax->title("Estimated Trajectory");
    data = data_history->GetData("rv_est");
  } else {
    std::cout << "State type must be either true or est." << std::endl;
    return;
  }

  // get x, y, z out of rv
  std::vector<double> x, y, z;
  for (auto& d : data) {
    x.push_back(d.GetData()(0));
    y.push_back(d.GetData()(1));
    z.push_back(d.GetData()(2));
  }
  ax->plot3(x, y, z, "b");

  // show figure
  show();
};

void PlotState(const std::shared_ptr<DataHistory> data_history,
               std::string state_type = "true") {
  using namespace matplot;

  if (state_type != "true" && state_type != "est") {
    std::cout << "State type must be either true or est." << std::endl;
    return;
  }

  std::vector<Timestamped<VectorXd>> data_rv;
  std::vector<Timestamped<VectorXd>> data_clk;

  if (state_type == "true") {
    data_rv = data_history->GetData("rv");
    data_clk = data_history->GetData("clk");
  } else {
    data_rv = data_history->GetData("rv_est");
    data_clk = data_history->GetData("clk_est");
  }

  // Plot the time history of each state, position in subplot 1, velocity in 2,
  // clock in 3
  auto fig = figure(true);
  auto ax1 = subplot(3, 1, 0);
  auto ax2 = subplot(3, 1, 1);
  auto ax3 = subplot(3, 1, 2);

  ax1->hold(on);
  ax1->grid(on);
  ax1->xlabel("Time [s]");
  ax1->ylabel("Position [km]");
  if (state_type == "true") {
    ax1->title("True Position");
  } else {
    ax1->title("Estimated Position");
  }

  ax2->hold(on);
  ax2->grid(on);
  ax2->xlabel("Time [s]");
  ax2->ylabel("Velocity [km/s]");
  if (state_type == "true") {
    ax2->title("True Velocity");
  } else {
    ax2->title("Estimated Velocity");
  }

  ax3->hold(on);
  ax3->grid(on);
  ax3->xlabel("Time [s]");
  ax3->ylabel("Clock [s]");
  if (state_type == "true") {
    ax3->title("True Clock");
  } else {
    ax3->title("Estimated Clock");
  }

  // get x, y, z out of rv
  std::vector<double> t, x, y, z, vx, vy, vz, clk_bias, clk_drift;
  for (auto& d : data_rv) {
    t.push_back(d.GetTimestamp());
    x.push_back(d.GetData()(0));
    y.push_back(d.GetData()(1));
    z.push_back(d.GetData()(2));
    vx.push_back(d.GetData()(3));
    vy.push_back(d.GetData()(4));
    vz.push_back(d.GetData()(5));
  }
  for (auto& d : data_clk) {
    clk_bias.push_back(d.GetData()(0));
    clk_drift.push_back(d.GetData()(1));
  }

  ax1->plot(t, x, "b");
  ax1->plot(t, y, "r");
  ax1->plot(t, z, "g");

  ax2->plot(t, vx, "b");
  ax2->plot(t, vy, "r");
  ax2->plot(t, vz, "g");

  ax3->plot(t, clk_bias, "b");
  ax3->plot(t, clk_drift, "r");

  ax1->legend({"x [km]", "y [km]", "z [km]"});
  ax2->legend({"vx [km/s]", "vy [km/s]", "vz [km/s]"});
  ax3->legend({"clk_bias [s]", "clk_drift [s/s]"});
  show();
};

int main() { return 0; };

}  // namespace lupnt