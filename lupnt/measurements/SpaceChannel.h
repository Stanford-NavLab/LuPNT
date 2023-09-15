/**
 * @file SpaceChannel.h
 * @author Stanford NAV LAB
 * @brief Base spacechannel class (Under devlopment )
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <lupnt/measurements/Occultation.h>

#include <memory>

namespace LPT {

class Transmission;
class ICommDevice;

class SpaceChannel {
 public:
  /**
   * @brief Calculates the transmitter-receiver link budget.
   *
   * This function calculates the full link budget between a transmitter sytem
   * and a receiver system. It includes the transmitter power, the antenna gain,
   * the path loss, atmospheric attenuation, the receiver antenna gain, the
   * receiver system noise temperature, receiver system noise figures and
   * losses, and received carrier to noise ratio. It uses a set of system gain
   * and loss constants with a set of distances and angles relative to the
   * transmitter and receiver antennas to describe the link budget geometry.
   *
   * Several methods of calculation can be used depending on the input
   * specification: 1) Calculate CN0 (and other outputs) from both the transmit
   * and receive antenna patterns, the antenna angles, and the LOS magnitude. 2)
   * Calculate CN0 (and other outputs) from both the transmit and receive
   * antenna patterns, the antenna angles, and the LOS magnitude, and the CN0
   * value. 3) Calculate CN0 (and other outputs) from the transmit antenna
   * pattern, the transmit antenna angles, the LOS magnitude, and the CN0 value.
   *
   * Either antenna pattern can be 1D, using only the antenna elevation angles,
   * or 2D, using the antenna azimuth and elevation angles. 1D patterns will
   * ignore any supplied azimuth angle data. Both types of antenna angles must
   * be supplied for a 2D antenna. The transmit antenna and receive antenna
   * patterns can be different dimensions (i.e. one can be 1D while the other is
   * 2D).
   *
   * If the receive pattern and receiver angles arguments are all empty then an
   * ideal omnidirectional receive antenna with zero dB gain is assumed.
   *
   * Transmit or receive gains, At or Ar, are set to -100 dB if the angles are
   * outside the prescribed pattern. (This is not applied to the receive gains
   * with the omni receive option.) All dependent outputs (CN0, RP, or AP, as
   * applicable) are biased as well.
   *
   * @param link
   *
   * @see linkbudget.m for the original MATLAB code and documentation from the
   * ODTBX project. The original MATLAB code is licensed under the BSD 3-Clause
   * License.
   *
   * @todo Check copyright and license information.
   *
   * The following are the link budget equations:
   * AP = P_sv + At + Ad + Ae;
   * RP = AP + Ar + As;
   * CN0 = RP - (10 * log10(Ts)) + 228.6 + Nf + L;
   *
   */
  void computeLinkBudget(std::shared_ptr<ICommDevice> &tx,
                         std::shared_ptr<ICommDevice> &rx, double t_tx,
                         double t_rx, Transmission &transmission);
};
}  // namespace LPT
