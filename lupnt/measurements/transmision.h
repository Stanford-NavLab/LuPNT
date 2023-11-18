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

#include <Eigen/Dense>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <string>

namespace ad = autodiff;

namespace lupnt {
struct Transmission {
  // Clock time [s]
  double t_tx;
  double t_rx;

  // frequency [Hz]
  double freq;
  std::string freq_label;

  // Clock offset from GNSS time [s]
  double dt_tx;
  double dt_rx;

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

  int ID_tx;
  Eigen::Vector3d r_tx, v_tx, r_rx, v_rx;
};
}  // namespace lupnt
