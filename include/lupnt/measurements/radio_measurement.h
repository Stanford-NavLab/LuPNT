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

#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {

enum Modulation {
  BPSK,   // Binary Phase Shift Keying
  QPSK,   // Quadrature Phase Shift Keying
  OQPSK,  // Offset Quadrature Phase Shift Keying
  GMSK,   // Gaussian Minimum Shift Keying
}

enum FrequencyBand {
  S,
  X,
  Ka
}

class RadioMeasurement {
 public:
  /**
   * @brief Compute the one-way range between two points
   *
   * @param r_tx   Transmitter position
   * @param r_rx   Receiver position
   * @param offset  Measurement offset
   * @return Real   One-way range
   */
  static Real ComputeOneWayRange(VecX r_tx, VecX r_rx, Real offset);

  /**
   * @brief
   *
   * @param r_tx
   * @param r_rx
   * @param dt_tx
   * @param dt_rx
   * @param offset
   * @return Real
   */
  static Real ComputePseudorange(VecX r_tx, VecX r_rx, Real dt_tx, Real dt_rx,
                                 Real offset);

  static Real ComputePseudorangerate(VecX r_tx, VecX r_rx, VecX v_tx, VecX v_rx,
                                     Real dt_tx_dot, Real dt_rx_dot,
                                     Real offset);

  static Real ComputeDopplerShift(VecX r_tx, VecX r_rx, VecX v_tx, VecX v_rx,
                                  Real dt_tx_dot, Real dt_rx_dot, Real f,
                                  Real offset);

  /**
   * @brief Compute the one-way range between two points considering light time
   * delay
   *   Reference: Grenfell MIT Ph.D. thesis, 2024  (A.2)
   *
   * @param rv_tx_tr   Transmitter position at reception time t_R (w.r.t to
   * central body)
   * @param rv_rx_tr   Receiver position at reception time t_R  (w.r.t to
   * central body)
   * @param dt_tx      Transmitter clock offset
   * @param dt_rx      Receiver clock offset
   * @param GM         Gravitational constant of the central body
   * @param offset    Measurement offset
   * @return Real      One-way pseudorange at t_R (clock offset error included)
   */
  static Real ComputeOneWayRangeLTR(VecX rv_tx_tr, VecX rv_rx_tr, Real dt_tx,
                                    Real dt_rx, double GM, Real hardware_delay);

  /**
   * @brief Compute the one-way range between two points considering light time
   * delay  (rx->tx->rx)
   *   Reference: Grenfell MIT Ph.D. thesis, 2024  (A.2)
   *
   * @param rv_target_tr   Transmitter position at reception time t_R (w.r.t to
   * central body)
   * @param rv_rx_tr   Receiver position at reception time t_R  (w.r.t to
   * central body)
   * @param GM         Gravitational constant of the central body
   * @param offset    Measurement offset
   * @return Real      Range at t_R
   */
  static Real ComputeTwoWayRangeLTR(VecX rv_target_tr, VecX rv_rx_tr, double GM,
                                    Real hardware_delay);

  /**
   * @brief Compute the PN regenerative range error for chip tracking loop
   * (CTL): suited for onboard processing
   *  Reference: "Pseudo-Noise (PN) Ranging Systems", Greenbook 2014
   *
   * @param PRC_N0   Carrier-to-noise ratio for the range clock (linear)
   * @param B_L      One-sided Chip tracking Loop noise bandwidth (usually
   * around 1Hz, 0.5Hz)
   * @param T_c       Chip period [s] = = 1/ (2 f_RC)
   * @return double range error [m]
   */
  static double ComputePnRangeErrorCTL(double PRC_N0, double B_L, double T_c);

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
  static double ComputePnRangeErrorOL(double PRC_N0, double T_I, double T_c);

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
  static double ComputeRangeRateErrorOneWay(
      double B_L_carrier, double f_C, double T_s, double T_I, double PT_N0,
      double sigma_y_1s, Modulation carrier_type = Modulation::BPSK);

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
  static double ComputeRangeRateErrorTwoWay(
      double B_L_carrier, double f_C, double T_s, double T_I, double PT_N0,
      double sigma_y_1s, double G, Modulation carrier_type = Modulation::BPSK);

  /**
   * @brief Get the (recommended) Transponder Turn Around Ratio for spacecraft
   *  // https://deepspace.jpl.nasa.gov/dsndocs/810-005/201/201B.pdf
   *
   * @param fbu  uplink frequency band
   * @param fbd  downlink frequency band
   * @return double  Turn around ratio
   */
  static double GetTransponderTurnAroundRatio(FrequencyBand fbu,
                                              FrequencyBand fbd);

  /**
   * @brief Compute the carrier loop signal-to-noise ratio
   *
   * @param PT_N0  downlink total signal power to noise spectral density ratio
   * [Hz]
   * @param B_L_carrier  carrier loop noise bandwidth [Hz]
   * @param T_s           period of the binary symbol [s]
   * @param carrier_type  modulation type
   *
   */
  static double ComputeCarrierLoopSNR(double PT_N0, double B_L_carrier,
                                      double T_s, Modulation carrier_type)
};

}  // namespace lupnt