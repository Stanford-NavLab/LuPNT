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

#pragma once

#include "lupnt/measurements/comm_device.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/body.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {

  /**
   * @brief Compute the one-way range between two points
   *
   * @param r_tx   Transmitter position
   * @param r_rx   Receiver position
   * @param offset  Measurement offset
   * @return Real   One-way range
   */
  Real ComputeOneWayRange(VecX r_tx, VecX r_rx, Real offset);

  /**
   * @brief Compute the pseudorange between two points
   *
   * @param r_tx  Transmitter position
   * @param r_rx  Receiver position
   * @param dt_tx  Transmitter clock offset
   * @param dt_rx  Receiver clock offset
   * @param offset  Measurement offset
   * @return Real  Pseudorange
   */
  Real ComputePseudorange(VecX r_tx, VecX r_rx, Real dt_tx, Real dt_rx, Real offset);

  /**
   * @brief Compute the pseudorangerate between two points
   *
   * @param r_tx  Transmitter position
   * @param r_rx  Receiver position
   * @param v_tx  Transmitter velocity
   * @param v_rx  Receiver velocity
   * @param dt_tx_dot  Transmitter clock offset rate
   * @param dt_rx_dot  Receiver clock offset rate
   * @param offset  Measurement offset
   * @return Real   Pseudorangerate
   */
  Real ComputePseudorangerate(VecX r_tx, VecX r_rx, VecX v_tx, VecX v_rx, Real dt_tx_dot,
                              Real dt_rx_dot, Real offset);

  /**
   * @brief Compute the Doppler shift between two points
   *
   * @param r_tx  Transmitter position
   * @param r_rx  Receiver position
   * @param v_tx  Transmitter velocity
   * @param v_rx  Receiver velocity
   * @param dt_tx_dot  Transmitter clock offset rate
   * @param dt_rx_dot  Receiver clock offset rate
   * @param offset  Measurement offset
   * @return Real  Doppler shift [Hz]
   */
  Real ComputeDopplerShift(VecX r_tx, VecX r_rx, VecX v_tx, VecX v_rx, Real dt_tx_dot,
                           Real dt_rx_dot, Real f, Real offset);

  /**
   * @brief Compute the one-way range between two points considering light time
   * delay
   *   Reference: Grenfell MIT Ph.D. thesis, 2024  (A.2)
   *
   * @param epoch_rx   Reception epoch (TAI) t_R
   * @param rv_tx  Transmitter position at reception time t_R (w.r.t to
   * central body)
   * @param rv_rx   Receiver position at reception time t_R  (w.r.t to
   * central body)
   * @param dt_tx      Transmitter clock offset
   * @param dt_rx      Receiver clock offset
   * @param tx_center_body  Central body of the transmitter
   * @param rx_center_body  Central body of the receiver
   * @param is_bodyfixed_tx  Flag to indicate if the transmitter state is
   * body-fixed
   * @param is_bodyfixed_rx  Flag to indicate if the receiver state is
   * body-fixed
   * @param hardware_delay  Hardware delay
   * @return Real      One-way pseudorange at t_R (clock offset error included)
   */
  Real ComputeOneWayRangeLTR(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx, Real dt_tx, Real dt_rx,
                             BodyData, BodyData rx_center_body, bool is_bodyfixed_tx,
                             bool is_bodyfixed_rx, Real hardware_delay);

  /**
   * @brief Compute the two-way range between two points considering light time
   * delay  (rx -> target -> rx)
   *   Reference: Grenfell MIT Ph.D. thesis, 2024  (A.2)
   *
   * @param epoch_rx   Reception epoch (TAI) t_R
   * @param rv_target_tr  Transmitter position at reception time t_R (w.r.t to
   * central body)
   * @param rv_rx_tr   Receiver position at reception time t_R  (w.r.t to
   * central body)
   * @param target_center_body  Central body of the target
   * @param rx_center_body  Central body of the receiver
   * @param is_bodyfixed_target  Flag to indicate if the target state is
   * body-fixed
   * @param is_bodyfixed_rx  Flag to indicate if the receiver state is
   * body-fixed
   * @param hardware_delay  Hardware delay
   * @return Real      One-way pseudorange at t_R (clock offset error included)
   */
  Real ComputeTwoWayRangeLTR(Real epoch_rx, Vec6 rv_target_tr, Vec6 rv_rx_tr,
                             BodyData target_center_body, BodyData rx_center_body,
                             bool is_bodyfixed_target, bool is_bodyfixed_rx, Real hardware_delay,
                             Real additional_delay = 0.0);

  /**
   * @brief Compute the one-way range rate between two points considering light
   * time delay
   *
   * @param epoch_rx Reception epoch (TAI) t_R
   * @param rv_tx_tr  Transmitter position at reception time t_R (w.r.t to
   * @param rv_rx_tr  Receiver position at reception time t_R (w.r.t to
   * @param dt_dot_tx  Transmitter clock offset rate
   * @param dt_dot_rx  Receiver clock offset rate
   * @param target_center_body  Central body of the target
   * @param rx_center_body  Central body of the receiver
   * @param is_bodyfixed_target  Flag to indicate if the target state is
   * body-fixed
   * @param is_bodyfixed_rx  Flag to indicate if the receiver state is
   * body-fixed
   * @param hardware_delay  Hardware delay [s]
   * @param T_I  Integration time [s]
   * @return Real
   */
  Real ComputeOneWayRangeRateLTR(Real epoch_rx, Vec6 rv_tx_tr, Vec6 rv_rx_tr, Real dt_dot_tx,
                                 Real dt_dot_rx, BodyData target_center_body,
                                 BodyData rx_center_body, bool is_bodyfixed_target,
                                 bool is_bodyfixed_rx, Real hardware_delay, double T_I);

  Real ComputeTwoWayRangeRateLTR(Real epoch_rx, Vec6 rv_target_tr, Vec6 rv_rx_tr,
                                 BodyData target_center_body, BodyData rx_center_body,
                                 bool is_bodyfixed_target, bool is_bodyfixed_rx,
                                 Real hardware_delay, double T_I);

  /**
   * @brief Compute the PN regenerative range error for chip tracking loop
   * (CTL): suited for onboard processing
   *  Reference: "Pseudo-Noise (PN) Ranging Systems", Greenbook 2014
   *
   * @param PRC_N0   Carrier-to-noise ratio for the range clock (linear)
   * @param B_L_CTL      One-sided Chip tracking Loop noise bandwidth
   * (usually around 1Hz, 0.5Hz)
   * @param T_c       Chip period [s] = = 1/ (2 f_RC)
   * @return double range error [m]
   */
  double ComputePnRangeErrorCTL(double PRC_N0, double B_L_CTL, double T_c,
                                Modulation modulation_type = Modulation::BPSK);

  /**
   *
   * @brief Compute the PN regenerative range error for open loop (OL) tracking
   * suited for ground stations
   * Reference: "Pseudo-Noise (PN) Ranging Systems", Greenbook 2014
   *
   * @param PRC_N0 Carrier-to-noise ratio for the range clock (linear)
   * @param T_I Integration time
   * @param T_c Chip period [s] = 1/ (2 f_RC)
   * @return double   range error [m]
   */
  double ComputePnRangeErrorOL(double PRC_N0, double T_I, double T_c,
                               Modulation modulation_type = Modulation::BPSK);

  /**
   * @brief Compute range rate error
   * Reference: https://deepspace.jpl.nasa.gov/dsndocs/810-005/202/202E.pdf
   *
   * @param B_L_carrier
   * @param f_C   downlink carrier frequency  [Hz]
   * @param T_s   period of the binary symbol [s]
   * @param T_I   integration time [s]
   * @param sigma_y_1s  Allan deviation of the receiver clock at 1s
   * @param PT_N0  [Hz]
   * @return double  range rate error [m/s]
   */
  double ComputeRangeRateErrorOneWay(double B_L_carrier, double f_C, double T_s, double T_I,
                                     double PT_N0, double sigma_y_1s,
                                     Modulation modulation_type = Modulation::BPSK,
                                     double m_R = 0.0);

  /**
   * @brief Compute range rate error for two-way ranging
   * Reference: https://deepspace.jpl.nasa.gov/dsndocs/810-005/202/202E.pdf
   *
   * @param B_L_carrier  carrier loop noise bandwidth [Hz]
   * @param f_C   downlink carrier frequency [Hz]
   * @param T_s   period of the binary symbol [s]
   * @param T_I   integration time [s]
   * @param sigma_y_1s  Allan deviation of the receiver clock at 1s
   * @param G     transponder turnaround ratio
   * @param PT_N0  downlink total signal power to noise spectral density ratio
   * [Hz]
   * @return double  range rate error [m/s]
   */
  double ComputeRangeRateErrorTwoWay(double B_L_carrier, double f_C, double T_s, double T_I,
                                     double PT_N0, double sigma_y_1s, double G,
                                     Modulation modulation_type = Modulation::BPSK,
                                     double m_R = 0.0);

}  // namespace lupnt
