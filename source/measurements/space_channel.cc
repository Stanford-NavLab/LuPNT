/**
 * @file SpaceChannel.cpp
 * @author Stanford NAV LAB
 * @brief Base spacechannel class (Under devlopment )
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/measurements/space_channel.h"

#include <iomanip>
#include <iostream>
#include <string>

#include "lupnt/agents/agent.h"
#include "lupnt/core/constants.h"
#include "lupnt/measurements/comm_device.h"

namespace lupnt {

void SpaceChannel::ComputeLinkBudget(std::shared_ptr<ICommDevice> &tx,
                                     std::shared_ptr<ICommDevice> &rx, double t,
                                     std::string time_fixed,
                                     Transmission &transmission) {
  // Transmitter and receiver positions and velocities
  double tau = 0.0;  // light time delay
  std::shared_ptr<CartesianOrbitState> rv_rx_gcrf;
  std::shared_ptr<CartesianOrbitState> rv_tx_gcrf;
  double t_rx = t;
  double t_tx = t;

  // Solve light time delay
  if (time_fixed == "rx") {
    tau = SolveLightTimeDelayRx(tx, rx, t);
    rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t);
    rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t - tau);
    t_rx = t;
    t_tx = t - tau;
  } else if (time_fixed == "tx") {
    tau = SolveLightTimeDelayTx(tx, rx, t);
    rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t - tau);
    rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t);
    t_rx = t + tau;
    t_tx = t;
  } else {
    // generate error
    std::cerr << "Error: Invalid time_fixed parameter" << std::endl;
  }

  // Commpute Occultations
}

double SpaceChannel::SolveLightTimeDelayRx(std::shared_ptr<ICommDevice> &tx,
                                           std::shared_ptr<ICommDevice> &rx,
                                           double t_rx) {
  // Transmitter and receiver positions and velocities
  double tau = 0.0;  // light time delay
  auto rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_rx);
  auto rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_rx - tau);

  // Compute Light time delay
  double tau_prev = 0.0;  // propagation time
  int n_iter = 0;
  int max_iter = 100;

  for (int n_iter = 0; n_iter < max_iter; n_iter++) {
    rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_rx - tau);
    double rho = (rv_tx_gcrf->r() - rv_rx_gcrf->r()).norm().val();
    tau = rho / C;
    if (fabs(tau - tau_prev) < 1e-12)
      break;
    else {
      tau_prev = tau;
    }
  }

  return tau;
}

double SpaceChannel::SolveLightTimeDelayTx(std::shared_ptr<ICommDevice> &tx,
                                           std::shared_ptr<ICommDevice> &rx,
                                           double t_tx) {
  // Transmitter and receiver positions and velocities
  double tau = 0.0;  // light time delay
  auto rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_tx + tau);
  auto rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_tx);

  // Compute Light time delay
  double tau_prev = 0.0;  // propagation time
  int n_iter = 0;
  int max_iter = 100;

  for (int n_iter = 0; n_iter < max_iter; n_iter++) {
    rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_tx + tau);
    double rho = (rv_tx_gcrf->r() - rv_rx_gcrf->r()).norm().val();
    tau = rho / C;
    if (fabs(tau - tau_prev) < 1e-12)
      break;
    else {
      tau_prev = tau;
    }
  }

  return tau;
}

double SpaceChannel::ComputeFreeSpaceLossdB(double dist, double lambda) {
  double path_loss = pow(4 * PI * dist / lambda, 2);
  double path_loss_dB = 10 * log10(path_loss);
  return path_loss_dB;
}

}  // namespace lupnt