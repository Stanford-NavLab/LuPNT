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

#include "lupnt/measurements/occultation.h"

namespace lupnt {

  class Transmission;
  class ICommDevice;

  class SpaceChannel {
  public:
    /**
     * @brief Calculates the transmitter-receiver link budget.
     *
     * The following are the link budget equations:
     * AP = P_sv + At + Ad + Ae;
     * RP = AP + Ar + As;
     * CN0 = RP - (10 * log10(Ts)) + 228.6 + Nf + L;
     *
     */
    void computeLinkBudget(std::shared_ptr<ICommDevice> &tx, std::shared_ptr<ICommDevice> &rx,
                           double t_tx, double t_rx, Transmission &transmission);
  };
}  // namespace lupnt
