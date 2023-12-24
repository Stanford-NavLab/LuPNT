/**
 * @file example_utils.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2023-12-23
 *
 * @copyright Copyright (c) 2023
 *
 */

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

FilterProcessNoiseFunction proc_noise_func = [](const VectorX& x, real t_curr,
                                                real t_end) {
  int clock_index = 6;
  double dt = (t_end - t_curr).val();
  double sigma_acc = 1e-13;

  MatrixXd Q = MatrixXd::Zero(state_size, state_size);

  Matrix6d Q_rv = Matrix6d::Zero();
  for (int i = 0; i < 3; i++) {
    Q_rv(i, i) = pow(dt, 3) / 3.0 * pow(sigma_acc, 2);
    Q_rv(i + 3, i + 3) = dt * pow(sigma_acc, 2);
    Q_rv(i, i + 3) = pow(dt, 2) / 2.0 * pow(sigma_acc, 2);
    Q_rv(i + 3, i) = pow(dt, 2) / 2.0 * pow(sigma_acc, 2);
  }
  Matrix2d Q_clk = GetClockProcessNoise(ClockModel::kMicrosemiCsac, dt);

  Q.block(0, 0, 6, 6) = Q_rv;
  Q.block(6, 6, 2, 2) = Q_clk;

  return Q;
};

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
  data_history->AddData("z_true", t, ekf->z_true);
  data_history->AddData("z_pred", t, ekf->z_pred);
  data_history->AddData("CN0", t, meas->GetCN0());

  data_history->AddData("vis_earth", t, meas->GetEarthOccultation());
  data_history->AddData("vis_moon", t, meas->GetMoonOccultation());
  data_history->AddData("vis_antenna", t, meas->GetMoonOccultation());
  data_history->AddData("vis_ionos", t, meas->GetMoonOccultation());

  // Moon spacecraft
  auto state = sat->GetCartesianGCRFStateAtEpoch(epoch);
  auto sate_mi = ConvertOrbitStateCoordSystem(state, epoch, CoordSystem::MI);
  auto state_gcrf =
      ConvertOrbitStateCoordSystem(state, epoch, CoordSystem::GCRF);
  data_history->AddData("rv_moon_mi", t, sate_mi->GetVector());
  data_history->AddData("rv_moon_gcrf", t, state_gcrf->GetVector());

  // Estimation
  data_history->AddData("rv", t, sat->GetOrbitState()->GetVector());
  data_history->AddData("rv_pred", t, ekf->xbar.head(6));
  data_history->AddData("rv_est", t, ekf->x.head(6));

  data_history->AddData("clk", t, sat->GetClockState().GetVector());
  data_history->AddData("clk_pred", t, ekf->xbar.tail(2));
  data_history->AddData("clk_est", t, ekf->x.tail(2));

  data_history->AddData("P_rv", t, ekf->P.diagonal().segment(0, 6));
  data_history->AddData("P_clk", t, ekf->P.diagonal().segment(6, 2));

  // GPS constellation
  for (int i = 0; i < gps_const->GetNumSatellites(); i++) {
    auto sate = gps_const->GetSatellite(i)->GetCartesianGCRFStateAtEpoch(epoch);
    auto sate_mi = ConvertOrbitStateCoordSystem(sate, epoch, CoordSystem::MI);
    auto state_gcrf =
        ConvertOrbitStateCoordSystem(sate, epoch, CoordSystem::GCRF);

    std::string name = "sat" + std::to_string(i);
    data_history->AddData(name + "_mi", t, sate->GetVector());
    data_history->AddData(name + "_gcrf", t, state_gcrf->GetVector());
  }

  // Bodies
  data_history->AddData(
      "earth_mi", t,
      CoordConverter::Convert(epoch, VectorX::Zero(6), CoordSystem::GCRF,
                              CoordSystem::MI));
  data_history->AddData(
      "moon_gcrf", t,
      CoordConverter::Convert(epoch, VectorX::Zero(6), CoordSystem::MI,
                              CoordSystem::GCRF));
};

void PrintProgress(double t, const std::shared_ptr<Spacecraft> sat, EKF* ekf) {
  auto x_est = ekf->x;
  auto x_true = sat->GetStateVector();

  double x_pos_err = (x_true.segment(0, 3) - x_est.segment(0, 3)).norm().val();
  double x_vel_err = (x_true.segment(3, 3) - x_est.segment(3, 3)).norm().val();
  double x_clk_bias_err = abs((x_true(6) - x_est(6)).val());

  std::cout.precision(3);
  std::cout << t / 3600 << " | " << x_pos_err << " | " << x_vel_err << " | "
            << x_clk_bias_err << std::endl;

  // true and estimated state
  std::cout << "x_true: " << std::endl << x_true << std::endl;
  std::cout << "x_est: " << std::endl << x_est << std::endl;
  std::cout << " " << std::endl;
};

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

}  // namespace lupnt