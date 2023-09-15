/**
 * @file RadioMeasurement.cpp
 * @author Stanford NAV LAB
 * @brief Class for Radionavigation measurements
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "RadioMeasurement.h"

namespace LPT {
ad::real RadioMeasurement::ComputePseudorange(ad::VectorXreal r_tx,
                                              ad::VectorXreal r_rx,
                                              ad::real dt_tx, ad::real dt_rx,
                                              ad::real offset) {
  // P_rx = rho_rx + c*(dt_rx(t_rx) - dt_tx(t_tx)) + I_rx + T_rx + eps_P
  ad::real rho_rx = (r_tx - r_rx).norm();
  ad::real P_rx = rho_rx + C * (dt_rx - dt_tx) + offset;
  return P_rx;
};

ad::real RadioMeasurement::ComputePseudorangerate(
    ad::VectorXreal r_tx, ad::VectorXreal r_rx, ad::VectorXreal v_tx,
    ad::VectorXreal v_rx, ad::real dt_tx_dot, ad::real dt_rx_dot,
    ad::real offset) {
  ad::VectorXreal e_rx = (r_tx - r_rx).normalized();
  ad::real prr = dot(e_rx, v_tx - v_rx) + C * (dt_rx_dot - dt_tx_dot) + offset;
  return prr;
};

ad::real RadioMeasurement::ComputeDopplerShift(
    ad::VectorXreal r_tx, ad::VectorXreal r_rx, ad::VectorXreal v_tx,
    ad::VectorXreal v_rx, ad::real dt_tx_dot, ad::real dt_rx_dot, ad::real f,
    ad::real offset) {
  ad::real f_D = -f / C *
                 ComputePseudorangerate(r_tx, r_rx, v_tx, v_rx, dt_tx_dot,
                                        dt_rx_dot, offset);
  return f_D;
};

}  // namespace LPT