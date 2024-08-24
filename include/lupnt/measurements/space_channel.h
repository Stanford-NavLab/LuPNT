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

  struct ITransmission;

  class ICommDevice;
  class Transmitter;
  class Receiver;

  class SpaceChannel {
  private:
    std::vector<NaifId> occult_bodies_;
    VecXd occult_alt_;

  public:
    SpaceChannel() = default;

    /**
     * @brief Construct a new Space Channel object
     *
     * @param occult_bodies  occulting bodies
     * @param occult_alt  occultation altitude
     */
    void SetOccultationBodies(std::vector<NaifId> occult_bodies, VecXd occult_alt) {
      occult_bodies_ = occult_bodies;
      occult_alt_ = occult_alt;
    }

    /**
     * @brief Compute the link budget
     *
     * @param txDevice  transmitter device
     * @param rxDevice  receiver device
     * @param t         epoch of transmission
     * @param time_fixed  time fixed at transmitter or receiver (tx or rx)
     * @param transmission  transmission object
     */
    ITransmission ComputeLinkBudget(std::shared_ptr<Transmitter> &txDevice,
                                    std::shared_ptr<Receiver> &rxDevice, Real t,
                                    std::string time_fixed);

    /**
     * @brief  Compute the link budget for a given data rate
     *
     * @param txDevice
     * @param rxDevice
     * @param t
     * @param time_fixed
     * @return ITransmission
     */
    ITransmission ComputeLinkBudgetDR(std::shared_ptr<Transmitter> &txDevice,
                                      std::shared_ptr<Receiver> &rxDevice, Real t,
                                      std::string time_fixed, double data_rate);

    /**
     * @brief Solve the light time delay at the receiver
     *
     * @param tx  transmitter device
     * @param rx  receiver device
     * @param t_rx  receiver time
     * @return Real light time delay
     */
    Real SolveLightTimeDelayRx(std::shared_ptr<Transmitter> &tx, std::shared_ptr<Receiver> &rx,
                               Real t_rx);

    /**
     * @brief Solve the light time delay at the transmitter
     *
     * @param tx   transmitter device
     * @param t_tx   transmitter time
     * @return double   light time delay
     */
    Real SolveLightTimeDelayTx(std::shared_ptr<Transmitter> &tx, std::shared_ptr<Receiver> &rx,
                               Real t_tx);

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
