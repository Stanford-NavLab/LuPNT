/**
 * @file gnss_measurement.cpp
 * @author Stanford NAV LAB
 * @brief Class that constructs GPS measurement data from channel observables
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "gnss_measurement.h"

#include <lupnt/numerics/math_utils.h>
#include <lupnt/physics/coord_converter.h>

#include "radio_measurement.h"

namespace lupnt {

GnssMeasurement::GnssMeasurement(const std::vector<Transmission> trans)
    : trans_store(trans),
      f(trans.size()),
      dt_tx(trans.size()),
      I_rx(trans.size()),
      T_rx(trans.size()),
      eps_P(trans.size()),
      r_tx(3, trans.size()),
      r_rx(3),
      CN0(trans.size()),
      vis_earth(trans.size()),
      vis_moon(trans.size()),
      vis_antenna(trans.size()),
      vis_atmos(trans.size()),
      vis_ionos(trans.size()),
      rho_rx(trans.size()),
      P_rx(trans.size()) {
  n_meas = trans.size();

  // Iterate over received transmissions and update dt_tx, I_rx, and T_rx
  int i = 0;
  for (auto tr : trans) {
    // Clock time
    t_rx = tr.t_rx;
    t_tx = tr.t_tx;

    // Clock offset
    dt_rx = tr.dt_rx;
    dt_tx[i] = tr.dt_tx;

    // Position
    r_rx = tr.r_rx;
    r_tx.col(i) = tr.r_tx;

    // Propagation
    I_rx[i] = tr.I_rx;
    T_rx[i] = tr.T_rx;
    CN0[i] = tr.CN0;
    eps_P[i] = 0.0;

    // Visibility
    vis_earth[i] = tr.vis_earth;
    vis_moon[i] = tr.vis_moon;
    vis_antenna[i] = tr.vis_antenna;
    vis_atmos[i] = tr.vis_atmos;
    vis_ionos[i] = tr.vis_ionos;

    f[i] = tr.freq;

    ID_tx.push_back(tr.ID_tx);

    i++;
  }
}

/**
 * @brief Create another measurement class instance with certain band
 *
 * @param freq_label
 * @return GnssMeasurement
 */
GnssMeasurement GnssMeasurement::ExtractSignal(std::string freq_label) {
  std::vector<Transmission> transmissions_freq;
  for (auto &tx : trans_store) {
    if (tx.freq_label == freq_label) {
      transmissions_freq.push_back(tx);
    }
  }

  return GnssMeasurement(transmissions_freq);
}

VectorX GnssMeasurement::ComputePseudorange(VectorX r_rx, real dt_rx) const {
  // P_rx = rho_rx + c*(dt_rx(t_rx) - dt_tx(t_tx)) + I_rx + T_rx + eps_P
  VectorX P_rx(r_tx.cols());
  for (int i = 0; i < r_tx.cols(); i++) {
    P_rx(i) = RadioMeasurement::ComputePseudorange(r_rx, r_tx.col(i), dt_tx[i],
                                                   dt_rx, I_rx(i) + T_rx(i));
  }
  return P_rx;
}

VectorX GnssMeasurement::GetPseudorange() {
  return ComputePseudorange(r_rx, dt_rx);
}

VectorX GnssMeasurement::GetPseudorange(double epoch, Vector6 rv_pred,
                                        Vector2 clk_pred, MatrixXd &H_pr) {
  auto func = [&, epoch, this](const Vector6 &rv_pred,
                               const Vector2 &clk_pred) {
    auto rv_pred_gcrf = CoordConverter::Convert(epoch, rv_pred, CoordSystem::MI,
                                                CoordSystem::GCRF);
    Vector3 r_rx = rv_pred_gcrf.head(3);
    real dt_rx = clk_pred(0);
    return ComputePseudorange(r_rx, dt_rx);
  };

  Vector6 rv_pred_ad = rv_pred.cast<double>();
  Vector2 clk_pred_ad = clk_pred.cast<double>();

  VectorX z_pr_pred(n_meas);
  H_pr = jacobian(func, wrt(rv_pred_ad, clk_pred_ad),
                  at(rv_pred_ad, clk_pred_ad), z_pr_pred);
  return z_pr_pred;
}

VectorX GnssMeasurement::GetCarrierPhase() {
  // phi_rx = c / lambda * (t_rx - t_tx) + c / lambda * (dt_rx(t_rx) -
  // dt_tx(t_tx)) + phi_rx_0 - phi_0 + N_rx + eps_phi

  VectorXd phi_rx =
      c * lambda * (t_rx - t_tx + dt_rx - dt_tx.array()).matrix() +
      (phi_rx_tx - phi_tx + N_rx + eps_phi);
  return phi_rx;
}

VectorX GnssMeasurement::GetPhaseRange() {
  // Phi_rx = c * (t_rx - t_tx) + c*(dt_rx(t_rx)  dt_tx(t_tx)) + lambda *
  // (phi_rx_0 - phi_0 + N_rx) + lambda*eps_Phi

  VectorX Phi_rx = c * (t_rx - t_tx + dt_rx - dt_tx.array()).matrix() +
                   lambda * (phi_rx_tx - phi_tx + N_rx + eps_Phi);
  return Phi_rx;
};

VectorX GnssMeasurement::GetDopplerShift() {
  // f_D = - f/c*((v_tx(t_tx) - v_rx(t_rx))^T * e_rx + c * dt_rx_dot(t_rx) -
  // c* dt_tx_dot(t_tx))) + eps_D

  VectorXd f_D = -f / c *
                     (((v_tx - v_rx) * e_rx).array() + c * dt_rx_dot -
                      c * dt_tx_dot.array())
                         .matrix() +
                 eps_D;
  return f_D;
}

VectorX GnssMeasurement::GetPseudorangeRate(const VectorX &r_rx_,
                                            const VectorX &v_rx_) {
  // f_D = - f/c*((v_tx(t_tx) - v_rx(t_rx))^T * e_rx + c * dt_rx_dot(t_rx) -
  // c* dt_tx_dot(t_tx))) + eps_D

  VectorXd f_D =
      (((v_tx - v_rx) * e_rx).array() + c * dt_rx_dot - c * dt_tx_dot.array())
          .matrix() +
      eps_D;
  return f_D;
}
}  // namespace lupnt