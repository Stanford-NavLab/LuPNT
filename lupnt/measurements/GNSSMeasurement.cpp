/**
 * @file GNSSMeasurement.cpp
 * @author Stanford NAV LAB
 * @brief Class that constructs GPS measurement data from channel observables
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "GNSSMeasurement.h"

#include <lupnt/numerics/MathUtils.h>
#include <lupnt/physics/CoordConverter.h>

#include "RadioMeasurement.h"

namespace LPT {

GNSSMeasurement::GNSSMeasurement(const std::vector<Transmission> trans)
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
 * @return GNSSMeasurement
 */
GNSSMeasurement GNSSMeasurement::ExtractSignal(std::string freq_label) {
  std::vector<Transmission> transmissions_freq;
  for (auto &tx : trans_store) {
    if (tx.freq_label == freq_label) {
      transmissions_freq.push_back(tx);
    }
  }

  return GNSSMeasurement(transmissions_freq);
}

ad::VectorXreal GNSSMeasurement::ComputePseudorange(ad::VectorXreal r_rx,
                                                    ad::real dt_rx) const {
  // P_rx = rho_rx + c*(dt_rx(t_rx) - dt_tx(t_tx)) + I_rx + T_rx + eps_P
  ad::VectorXreal P_rx(r_tx.cols());
  for (int i = 0; i < r_tx.cols(); i++) {
    P_rx(i) = RadioMeasurement::ComputePseudorange(r_rx, r_tx.col(i), dt_tx[i],
                                                   dt_rx, I_rx(i) + T_rx(i));
  }
  return P_rx;
}

ad::VectorXreal GNSSMeasurement::GetPseudorange() {
  return ComputePseudorange(r_rx, dt_rx);
}

ad::VectorXreal GNSSMeasurement::GetPseudorange(double epoch,
                                                ad::Vector6real rv_pred,
                                                ad::Vector2real clk_pred,
                                                Eigen::MatrixXd &H_pr) {
  auto func = [&, epoch, this](const ad::Vector6real &rv_pred,
                               const ad::Vector2real &clk_pred) {
    auto rv_pred_gcrf = CoordConverter::Convert(rv_pred, epoch, CoordSystem::MI,
                                                CoordSystem::GCRF);
    ad::Vector3real r_rx = rv_pred_gcrf.head(3);
    ad::real dt_rx = clk_pred(0);
    return ComputePseudorange(r_rx, dt_rx);
  };

  ad::Vector6real rv_pred_ad = toEigen(rv_pred);
  ad::Vector2real clk_pred_ad = toEigen(clk_pred);

  ad::VectorXreal z_pr_pred(n_meas);
  H_pr = ad::jacobian(func, ad::wrt(rv_pred_ad, clk_pred_ad),
                      ad::at(rv_pred_ad, clk_pred_ad), z_pr_pred);
  return z_pr_pred;
}

ad::VectorXreal GNSSMeasurement::GetCarrierPhase() {
  // phi_rx = c / lambda * (t_rx - t_tx) + c / lambda * (dt_rx(t_rx) -
  // dt_tx(t_tx)) + phi_rx_0 - phi_0 + N_rx + eps_phi

  Eigen::VectorXd phi_rx =
      c * lambda * (t_rx - t_tx + dt_rx - dt_tx.array()).matrix() +
      (phi_rx_tx - phi_tx + N_rx + eps_phi);
  return phi_rx;
}

ad::VectorXreal GNSSMeasurement::GetPhaseRange() {
  // Phi_rx = c * (t_rx - t_tx) + c*(dt_rx(t_rx)  dt_tx(t_tx)) + lambda *
  // (phi_rx_0 - phi_0 + N_rx) + lambda*eps_Phi

  ad::VectorXreal Phi_rx = c * (t_rx - t_tx + dt_rx - dt_tx.array()).matrix() +
                           lambda * (phi_rx_tx - phi_tx + N_rx + eps_Phi);
  return Phi_rx;
};

ad::VectorXreal GNSSMeasurement::GetDopplerShift() {
  // f_D = - f/c*((v_tx(t_tx) - v_rx(t_rx))^T * e_rx + c * dt_rx_dot(t_rx) -
  // c* dt_tx_dot(t_tx))) + eps_D

  Eigen::VectorXd f_D = -f / c *
                            (((v_tx - v_rx) * e_rx).array() + c * dt_rx_dot -
                             c * dt_tx_dot.array())
                                .matrix() +
                        eps_D;
  return f_D;
}

ad::VectorXreal GNSSMeasurement::GetPseudorangeRate(
    const ad::VectorXreal &r_rx_, const ad::VectorXreal &v_rx_) {
  // f_D = - f/c*((v_tx(t_tx) - v_rx(t_rx))^T * e_rx + c * dt_rx_dot(t_rx) -
  // c* dt_tx_dot(t_tx))) + eps_D

  Eigen::VectorXd f_D =
      (((v_tx - v_rx) * e_rx).array() + c * dt_rx_dot - c * dt_tx_dot.array())
          .matrix() +
      eps_D;
  return f_D;
}
}  // namespace LPT