/**
 * @file gnssChannel.cpp
 * @author Stanford NAV LAB
 * @brief Gnss Channel
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "lupnt/measurements/gnss_channel.h"

#include "lupnt/agents/agent.h"
#include "lupnt/core/constants.h"
#include "lupnt/measurements/gnss_receiver.h"
#include "lupnt/measurements/gnss_transmitter.h"
#include "lupnt/physics/orbit_state.h"
#include "lupnt/physics/spice_interface.h"

#define DEBUG_TRANSMISSIONS 0

namespace lupnt {
  /**
   * @brief Receiver receives all available Gnss signals
   *
   * @param rx
   * @param t
   * @return std::vector<GnssTransmission>
   */
  std::vector<GnssTransmission> GnssChannel::Receive(GnssReceiver &rx, double t) {
    std::vector<GnssTransmission> received_transs;  // create an empty vector

    // Messages from other comms systems that can generate Gnss messages
    for (auto &tx : tx_devices) {
      // Solve light time delay
      double tau = 0.0;  // light time delay
      CartesianOrbitState rv_rx_gcrf = rx.GetAgent()->GetCartesianGCRFStateAtEpoch(t);
      CartesianOrbitState rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t - tau);

      // Compute Light time delay
      double tau_prev = 0.0;  // propagation time
      int n_iter = 0;
      int max_iter = 100;

      for (int n_iter = 0; n_iter < max_iter; n_iter++) {
        rv_tx_gcrf = tx->GetAgent()->GetCartesianGCRFStateAtEpoch(t - tau);
        double rho = (rv_tx_gcrf.r() - rv_rx_gcrf.r()).norm().val();
        tau = rho / C;
        if (fabs(tau - tau_prev) < 1e-12)
          break;
        else {
          tau_prev = tau;
        }
      }

      // Transmission and reception times
      double t_rx = t;
      double t_tx = t - tau;

      // Convert to Moon Inertial frame
      auto rv_tx_mi = ConvertOrbitStateFrame(rv_tx_gcrf, t_tx, Frame::MOON_CI);
      auto rv_rx_mi = ConvertOrbitStateFrame(rv_rx_gcrf, t_rx, Frame::MOON_CI);

      // Occultation
      std::string tx_planet = "";
      std::map<std::string, bool> vis = Occultation::ComputeOccultationGnss(
          rv_tx_gcrf.r().cast<double>(), rv_tx_mi.r().cast<double>(), rv_rx_gcrf.r().cast<double>(),
          rv_rx_mi.r().cast<double>(), tx_planet);

      if (vis["EARTH"] || vis["MOON"]) continue;  // quit if occulted

      // Transmitter and Receiver Antenna gain
      double At = tx->GetTransmitterAntennaGain(t_tx, rv_tx_gcrf.r().cast<double>(),
                                                rv_rx_gcrf.r().cast<double>());
      double Ar = rx.GetReceiverAntennaGain(t_rx, rv_tx_gcrf.r().cast<double>(),
                                            rv_rx_gcrf.r().cast<double>());

      // Generate transmission
      GnssTransmission trans = tx->GenerateTransmission(t_tx);
      double d = (rv_tx_gcrf.r() - rv_rx_gcrf.r()).norm().val();

      // Link budget
      for (size_t freq_idx = 0; freq_idx < tx->freq_list.size(); freq_idx++) {
        std::string freq_name = tx->freq_list[freq_idx];
        double freq = tx->freq_map[freq_name];
        double Ad = 20.0 * log10((C / freq) / (4.0 * PI * d));
        double scalars = tx->P_tx + rx.rx_param_.Ae + rx.rx_param_.As
                         - (10.0 * log10(rx.rx_param_.Tsys)) + 228.6 + rx.rx_param_.L;
        trans.CN0 = At + Ar + Ad + scalars;

        if (At <= -499.0 || vis["earth"] || vis["moon"] || trans.CN0 < rx.rx_param_.CN0threshold) {
          // not visible
          continue;
          trans.CN0 = NAN;
          trans.AP = NAN;
          trans.RP = NAN;
          trans.vis_antenna = true;
        } else {
          trans.AP = tx->P_tx + At + Ad + rx.rx_param_.Ae;
          trans.RP = trans.AP + Ar + rx.rx_param_.As;
          trans.vis_antenna = false;
        }

        // TX
        trans.t_tx = t_tx;
        trans.freq = freq;
        trans.freq_label = freq_name;
        trans.chip_rate = tx->rc_map[freq_name];
        trans.dt_tx = 0.0;      // Todo: Get this from ephemeris
        trans.dt_tx_dot = 0.0;  // Todo: Get this from ephemeris
        trans.r_tx = rv_tx_gcrf.r().cast<double>();
        trans.v_tx = rv_tx_gcrf.v().cast<double>();

        // Channel
        trans.I_rx = 0.0;
        trans.T_rx = 0.0;
        trans.vis_atmos = vis["atmos"];
        trans.vis_ionos = vis["ionos"];
        trans.vis_earth = vis["earth"];
        trans.vis_moon = vis["moon"];

        trans.ID_tx = tx->GetPRN();

        // RX
        trans.t_rx = t_rx;
        trans.dt_rx = rx.GetAgent()->GetClockState().GetValue(0).val();
        trans.dt_rx_dot = rx.GetAgent()->GetClockState().GetValue(1).val();
        trans.r_rx = rv_rx_gcrf.r().cast<double>();
        trans.v_rx = rv_rx_gcrf.v().cast<double>();

        // receiver chip param
        trans.gnssr_param = rx.gnssr_param_;

        received_transs.push_back(trans);
      }
    }
    return received_transs;
  }  // Receive

}  // namespace lupnt
