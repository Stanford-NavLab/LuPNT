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
#include "lupnt/measurements/transmission.h"
#include "lupnt/physics/occultation.h"

namespace lupnt {

  ITransmission SpaceChannel::ComputeLinkBudget(std::shared_ptr<Transmitter> &tx,
                                                std::shared_ptr<Receiver> &rx, Real t,
                                                std::string time_fixed) {
    ITransmission trans;  // create an empty vector

    // Register the transmitter and receiver
    trans.tx = tx;
    trans.rx = rx;

    // Transmitter and receiver positions and velocities
    Real tau = 0.0;  // light time delay
    CartesianOrbitState rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t);
    CartesianOrbitState rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t);
    Real t_rx = t;
    Real t_tx = t;

    // Solve light time delay
    if (time_fixed == "rx") {
      tau = SolveLightTimeDelayRx(tx, rx, t);
      rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t - tau);
      t_rx = t;
      t_tx = t - tau;
    } else if (time_fixed == "tx") {
      tau = SolveLightTimeDelayTx(tx, rx, t);
      rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t - tau);
      t_rx = t + tau;
      t_tx = t;
    } else {
      // generate error
      std::cerr << "Error: Invalid time_fixed parameter" << std::endl;
    }

    // Commpute Occultations
    bool vis_all = true;
    std::map<std::string, bool> vis_occult;

    if (occult_bodies_.size() > 0) {
      Real epoch = (t_tx + t_rx) / 2.0;
      vis_occult
          = Occultation::ComputeOccultation(epoch, rv_tx_gcrf.r(), rv_rx_gcrf.r(), Frame::GCRF,
                                            Frame::GCRF, occult_bodies_, occult_alt_);
      vis_all = vis_occult["all"];
    }

    // Link Budget
    double At = tx->GetTransmitterAntennaGain(t_tx.val(), rv_tx_gcrf.r().cast<double>(),
                                              rv_rx_gcrf.r().cast<double>());
    double Ar = rx->GetReceiverAntennaGain(t_rx.val(), rv_tx_gcrf.r().cast<double>(),
                                           rv_rx_gcrf.r().cast<double>());

    double dist = (rv_tx_gcrf.r() - rv_rx_gcrf.r()).norm().val();
    double lambda = C / tx->freq_tx;
    double fsl_loss_dB = ComputeFreeSpaceLossdB(dist, lambda);

    double EIRP_dB = tx->P_tx + At;
    double G_T_rx_dB = Ar - 10.0 * log10(rx->rx_param_.Tsys);
    double loss = rx->rx_param_.Ae + rx->rx_param_.As + rx->rx_param_.L;  // sum of lossess (minus)
    double CN0 = EIRP_dB - fsl_loss_dB + 228.6 + G_T_rx_dB + loss;
    double CN = CN0 - 10.0 * log10(tx->bandwidth);
    double RP = CN0 - 228.6 + 10.0 * log10(rx->rx_param_.Tsys);  // Received Power
    double RP_N0 = RP - 10.0 * log10(rx->rx_param_.Tsys);        // Received Power Noise

    // register the values
    trans.EIRP = EIRP_dB;
    trans.G_T = G_T_rx_dB;
    trans.CN0 = CN0;
    trans.CN0_linear = pow(10, CN0 / 10.0);
    trans.vis_occult = vis_occult;
    trans.vis_all = vis_all;

    return trans;
  }

  Real SpaceChannel::SolveLightTimeDelayRx(std::shared_ptr<Transmitter> &tx,
                                           std::shared_ptr<Receiver> &rx, Real t_rx) {
    // Transmitter and receiver positions and velocities
    Real tau = 0.0;  // light time delay
    auto rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_rx);
    auto rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_rx - tau);

    // Compute Light time delay
    Real tau_prev = 0.0;  // propagation time
    int max_iter = 100;
    double rho = 0.0;

    for (int n_iter = 0; n_iter < max_iter; n_iter++) {
      rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_rx - tau);
      rho = (rv_tx_gcrf.r() - rv_rx_gcrf.r()).norm().val();
      tau = rho / C;
      if (fabs(tau.val() - tau_prev.val()) < 1e-12)
        break;
      else {
        tau_prev = tau;
      }
    }

    return tau;
  }

  Real SpaceChannel::SolveLightTimeDelayTx(std::shared_ptr<Transmitter> &tx,
                                           std::shared_ptr<Receiver> &rx, Real t_tx) {
    // Transmitter and receiver positions and velocities
    Real tau = 0.0;  // light time delay
    auto rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_tx + tau);
    auto rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_tx);

    // Compute Light time delay
    Real tau_prev = 0.0;  // propagation time
    int n_iter = 0;
    int max_iter = 100;
    double rho = 0.0;

    for (int n_iter = 0; n_iter < max_iter; n_iter++) {
      rv_rx_gcrf = rx->GetAgent()->GetCartesianGCRFStateAtEpoch(t_tx + tau);
      rho = (rv_tx_gcrf.r() - rv_rx_gcrf.r()).norm().val();
      tau = rho / C;
      if (fabs(tau.val() - tau_prev.val()) < 1e-12)
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
