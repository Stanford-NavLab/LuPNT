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

#include "lupnt/measurements/radio_measurement.h"

namespace lupnt {
  Real RadioMeasurement::ComputePseudorange(VecX r_tx, VecX r_rx, Real dt_tx, Real dt_rx,
                                            Real offset) {
    // P_rx = rho_rx + c*(dt_rx(t_rx) - dt_tx(t_tx)) + I_rx + T_rx + eps_P
    Real rho_rx = (r_tx - r_rx).norm();
    Real P_rx = rho_rx + C * (dt_rx - dt_tx) + offset;
    return P_rx;
  };

  Real RadioMeasurement::ComputePseudorangerate(VecX r_tx, VecX r_rx, VecX v_tx, VecX v_rx,
                                                Real dt_tx_dot, Real dt_rx_dot, Real offset) {
    VecX e_rx = (r_tx - r_rx).normalized();
    Real prr = e_rx.dot(v_tx - v_rx) + C * (dt_rx_dot - dt_tx_dot) + offset;
    return prr;
  };

  Real RadioMeasurement::ComputeDopplerShift(VecX r_tx, VecX r_rx, VecX v_tx, VecX v_rx,
                                             Real dt_tx_dot, Real dt_rx_dot, Real f, Real offset) {
    Real f_D
        = -f / C * ComputePseudorangerate(r_tx, r_rx, v_tx, v_rx, dt_tx_dot, dt_rx_dot, offset);
    return f_D;
  };

}  // namespace lupnt
