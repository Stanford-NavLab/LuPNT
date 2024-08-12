/**
 * @file space_channel.h
 * @author Stanford NAV LAB
 * @brief Base spacechannel class (Under devlopment )
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <memory>

#include "lupnt/physics/occultation.h"

namespace lupnt {

class Transmission;
class ICommDevice;

class SpaceChannel {
 public:
  /**
   * @brief Compute the link budget
   *
   * @param txDevice  transmitter device
   * @param rxDevice  receiver device
   * @param t         time
   * @param time_fixed  time fixed at transmitter or receiver (tx or rx)
   * @param transmission  transmission object
   */
  void ComputeLinkBudget(std::shared_ptr<ICommDevice> &txDevice,
                         std::shared_ptr<ICommDevice> &rxDevice, double t,
                         std::string time_fixed, Transmission &transmission);

  /**
   * @brief Solve the light time delay at the receiver
   *
   * @param tx  transmitter device
   * @param rx  receiver device
   * @param t_rx  receiver time
   * @return double  light time delay
   */
  double SolveLightTimeDelayRx(std::shared_ptr<ICommDevice> &tx,
                               std::shared_ptr<ICommDevice> &rx, double t_rx);

  /**
   * @brief Solve the light time delay at the transmitter
   *
   * @param tx   transmitter device
   * @param t_tx   transmitter time
   * @return double   light time delay
   */
  double SolveLightTimeDelayTx(std::shared_ptr<ICommDevice> &tx,
                               std::shared_ptr<ICommDevice> &rx, double t_tx);

  /**
   * @brief Compute the free space loss
   *
   * @param dist Distance between the transmitter and receiver [km]
   * @param lambda Wavelength of the signal [km]
   * @return double Free space loss [dB]
   */
  double ComputeFreeSpaceLossdB(double dist, double lambda);
};

}  // namespace lupnt
