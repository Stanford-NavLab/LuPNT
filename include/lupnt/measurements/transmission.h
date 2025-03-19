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

#include <map>
#include <string>
#include <vector>

#include "lupnt/measurements/comm_device.h"
#include "lupnt/measurements/gnss_receiver_param.h"

namespace lupnt {

  struct ITransmission {
    // Clock time [s]
    double t_tx;
    double t_rx;

    // frequency [Hz]
    double freq;
    std::string freq_label;

    // Position and velocity
    Vec3d r_tx, v_tx;
    Vec3d r_rx, v_rx;

    // Clock offset [s]
    double dt_tx, dt_rx;
    double dt_tx_dot, dt_rx_dot;

    // link budget
    double EIRP;        // Equivalent isotropic radiated power [dBW]
    double G_T;         // Transmit antenna gain / Noise temperature [dB/K]
    double CN0;         // Carrier-to-noise density [dB-Hz]
    double CN0_linear;  // Carrier-to-noise density [dB-Hz]

    // Signal power [W]
    double AP;
    double RP;

    // Transmitter and receiver
    std::shared_ptr<Transmitter> tx;
    std::shared_ptr<Receiver> rx;

    // Agents
    bool is_tx_gs;
    bool is_rx_gs;
    bool is_tx_bodyfixed;
    bool is_rx_bodyfixed;

    // visibility
    bool vis_all;
    std::map<std::string, bool> vis_occult;

    int ID_tx;
  };

  struct GnssTransmission : ITransmission {
    // TX
    int ID_tx;

    // Channel
    double I_rx;  // ionospheric delay [s]
    double T_rx;  // tropospheric delay [s]
    bool vis_atmos;
    bool vis_ionos;
    bool vis_earth;
    bool vis_moon;
    bool vis_antenna;

    // receiver chip param
    double chip_rate;  // receiver chip rate [Hz]
    GnssReceiverParam gnssr_param;
  };

}  // namespace lupnt
