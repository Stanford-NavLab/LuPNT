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

#include "radio_measurement.h"

namespace lupnt {
real RadioMeasurement::ComputePseudorange(VectorX r_tx, VectorX r_rx,
                                          real dt_tx, real dt_rx, real offset) {
  // P_rx = rho_rx + c*(dt_rx(t_rx) - dt_tx(t_tx)) + I_rx + T_rx + eps_P
  real rho_rx = (r_tx - r_rx).norm();
  real P_rx = rho_rx + C * (dt_rx - dt_tx) + offset;
  return P_rx;
};

real RadioMeasurement::ComputePseudorangerate(VectorX r_tx, VectorX r_rx,
                                              VectorX v_tx, VectorX v_rx,
                                              real dt_tx_dot, real dt_rx_dot,
                                              real offset) {
  VectorX e_rx = (r_tx - r_rx).normalized();
  real prr = e_rx.dot(v_tx - v_rx) + C * (dt_rx_dot - dt_tx_dot) + offset;
  return prr;
};

real RadioMeasurement::ComputeDopplerShift(VectorX r_tx, VectorX r_rx,
                                           VectorX v_tx, VectorX v_rx,
                                           real dt_tx_dot, real dt_rx_dot,
                                           real f, real offset) {
  real f_D = -f / C *
             ComputePseudorangerate(r_tx, r_rx, v_tx, v_rx, dt_tx_dot,
                                    dt_rx_dot, offset);
  return f_D;
};

}  // namespace lupnt