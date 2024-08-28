/**
 * @file comm_utils.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-08-21
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

namespace lupnt {

  enum Modulation {
    Residual,  // Residual carrier
    // Suppressed carrier
    BPSK,     // Binary Phase Shift Keying
    QPSK,     // Quadrature Phase Shift Keying
    OQPSK,    // Offset Quadrature Phase Shift Keying
    GMSK,     // Gaussian Minimum Shift Keying
    GMSK_PN,  // GMSK with PN modulation
  };

  enum FrequencyBand { UHF, L, S, Cband, X, Ku, K, Ka };

  /**
   * @brief Compute the Approximate BER for different Carrier modulation
   *  https://www.unilim.fr/pages_perso/vahid/notes/ber_awgn.pdf
   * @param PT_N0
   * @param Modulation
   * @return double
   */
  double ComputeBER(double EbN0, Modulation modulation_type);

  /**
   * @brief Compute Es/N0 from Eb/N0
   *
   * @param EbN0  Energy per bit to noise spectral density ratio
   * @param modulation_order  Modulation order (QPSK: 2, BPSK: 1)
   * @param coding_rate  Coding rate
   *
   */
  double ComputeEsN0(double EbN0, double modulation_order, double coding_rate);

  /**
   * @brief Compute the number of bits per symbol
   *
   * @param modulation_type  Modulation type
   * @return double
   */
  double BitsPerSymbol(Modulation modulation_type);
  /**
   * @brief Get the Frequency Band object
   *
   * @param f   frequency [Hz]
   * @return FrequencyBand
   */
  FrequencyBand GetFrequencyBand(double f);

  /**
   * @brief Get the (recommended) Transponder Turn Around Ratio for spacecraft
   *  // https://deepspace.jpl.nasa.gov/dsndocs/810-005/201/201B.pdf
   *
   * @param fbu  uplink frequency band
   * @param fbd  downlink frequency band
   * @return double  Turn around ratio
   */
  double GetTransponderTurnAroundRatio(FrequencyBand fbu, FrequencyBand fbd);

  /**
   * @brief Compute the carrier loop signal-to-noise ratio
   *
   * @param PT_N0  downlink total signal power to noise spectral density ratio
   * [Hz]
   * @param B_L_carrier  carrier loop noise bandwidth [Hz]
   * @param T_s           period of the binary symbol [s]
   * @param modulation_type  modulation type
   * @param m_R           modulation index (only for residual carrier)
   *
   */
  double ComputeCarrierLoopSNR(double PT_N0, double B_L_carrier, double T_s,
                               Modulation modulation_type, double m_R);

}  // namespace lupnt
