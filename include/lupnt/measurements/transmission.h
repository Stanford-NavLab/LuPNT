/**
 * @file transmission.h
 * @author Stanford NAV LAB
 * @brief Signal transmission data
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <string>

#include "gnss_receiver_param.h"

namespace lupnt {

  struct Transmission {
    // Clock time [s]
    double t_tx;
    double t_rx;

    // frequency [Hz]
    double freq;
    std::string freq_label;

    // Clock offset from Gnss time [s]
    double dt_tx;
    double dt_rx;

    // Clock drift from Gnss time [s/s]
    double dt_tx_dot;
    double dt_rx_dot;

    double I_rx;  // Ionospheric delay [m] (n_satellites * n_bands)
    double T_rx;  // Tropospheric delay [m] (n_satellites)
    double CN0;   // Carrier-to-noise density [dB-Hz] (n_satellites * n_bands)

    double AP;
    double RP;

    bool vis_earth;
    bool vis_moon;
    bool vis_antenna;
    bool vis_atmos;
    bool vis_ionos;

    GnssReceiverParam gnssr_param;

    double chip_rate;  // receiver chip rate [Hz]

    int ID_tx;
    Vec3d r_tx, v_tx, r_rx, v_rx;
  };
}  // namespace lupnt
