/**
 * @file link_measurement.cc
 * @author Stanford Nav Lab
 * @brief
 * @version 0.1
 * @date 2024-08-20
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "lupnt/measurements/link_measurement.h"

#include "lupnt/core/constants.h"
#include "lupnt/measurements/radio_measurement.h"
#include "lupnt/measurements/space_channel.h"
#include "lupnt/physics/body.h"

namespace lupnt {

  LinkMeasurement::LinkMeasurement(std::vector<NaifId> occult_bodies, VecXd occult_alt,
                                   Real hardware_delay)
      : hardware_delay_(hardware_delay), occult_bodies_(occult_bodies), occult_alt_(occult_alt) {}

  void LinkMeasurement::SetLinkParams() {
    ITransmission trans, trans_u;
    double G = 1.0;

    if (one_way_generated_) {
      trans = trans_ow_;

    } else if (two_way_generated_) {
      trans = trans_tw_[1];    // Downlink
      trans_u = trans_tw_[0];  // Uplink
    }

    // Link Parameters
    BodyData tx_center_body = GetBodyData(trans.tx->GetAgent()->GetBodyId());
    BodyData rx_center_body = GetBodyData(trans.rx->GetAgent()->GetBodyId());
    linkparams_.tx_center_body = tx_center_body;
    linkparams_.rx_center_body = rx_center_body;
    linkparams_.is_bodyfixed_tx = trans.tx->GetAgent()->IsBodyFixed();
    linkparams_.is_bodyfixed_rx = trans.rx->GetAgent()->IsBodyFixed();

    // Signal parameters
    ReceiverParam rx_param = trans.rx->rx_param_;
    linkparams_.freq = trans.tx->freq_tx;
    linkparams_.B_L_chip = rx_param.B_L_chip;
    linkparams_.B_L_carrier = rx_param.B_L_carrier;
    linkparams_.Tc = rx_param.modulation_type;
    linkparams_.T_I_doppler = rx_param.T_I_doppler;
    linkparams_.T_I_range = rx_param.T_I_range;

    // Singals
    linkparams_.CN0_linear = trans.CN0_linear;
  }

  /********************** One way Link ***************************/

  void LinkMeasurement::GenerateOneWayLink(Real epoch, std::shared_ptr<Transmitter> &tx,
                                           std::shared_ptr<Receiver> &rx, std::string txrx) {
    SpaceChannel sc = SpaceChannel();
    sc.SetOccultationBodies(occult_bodies_, occult_alt_);
    Real t_tx_d, t_rx_d;

    Real delay_tx = hardware_delay_;
    Real delay_rx = hardware_delay_;

    if (txrx == "rx") {
      t_rx_d = epoch;
      trans_ow_ = sc.ComputeLinkBudget(tx, rx, t_rx_d - delay_rx, "rx");
      t_tx_d = trans_ow_.t_tx - delay_tx;
    } else if (txrx == "tx") {
      t_tx_d = epoch;
      trans_ow_ = sc.ComputeLinkBudget(tx, rx, t_tx_d + delay_tx, "tx");
      t_rx_d = trans_ow_.t_rx + delay_rx;
    }
    one_way_generated_ = true;

    // Set time
    epoch_tx_true_ = t_tx_d;
    epoch_rx_true_ = t_rx_d;

    // Register the agents
    agents_.push_back(tx->GetAgent());
    agents_.push_back(rx->GetAgent());

    // Link Parameters
    SetLinkParams();
  }

  VecX LinkMeasurement::GetTrueOneWayLinkMeasurement(std::vector<LinkMeasurementType> meas_types) {
    if (!one_way_generated_) {
      std::cerr << "Error: One way link not generated" << std::endl;
      return VecX::Zero(0);
    }

    VecX rv_tx = agents_[0]->GetRvStateAtEpoch(epoch_rx_true_);
    VecX rv_rx = agents_[0]->GetRvStateAtEpoch(epoch_rx_true_);
    VecX clk_tx = agents_[1]->GetClockStateVecAtEpoch(epoch_rx_true_);
    VecX clk_rx = agents_[1]->GetClockStateVecAtEpoch(epoch_rx_true_);

    MatXd H_ow_rx(2, 8);  // temporary, won't be used

    return GetOneWayLinkMeasurement(epoch_rx_true_, rv_tx, rv_rx, clk_tx, clk_rx, H_ow_rx,
                                    hardware_delay_, meas_types, true, false);
  }

  VecX LinkMeasurement::GetOneWayLinkMeasurement(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx,
                                                 Vec2 clk_rx, MatXd H_ow_rx, Real hardware_delay,
                                                 std::vector<LinkMeasurementType> meas_types,
                                                 bool with_noise, bool with_jacobian) {
    Real rho_ow, rho_ow_rate;
    double sigma_ow = 0.0;
    double sigma_ow_rate = 0.0;

    int meas_size = meas_types.size();
    int state_size = state_size_ow_;

    VecX z(meas_size);
    H_ow_rx.resize(meas_size, state_size);
    MatXd H_ow_range = MatXd::Zero(1, state_size);
    MatXd H_ow_rangerate = MatXd::Zero(1, state_size);

    int idx = 0;

    if (meas_size == 0) {
      return VecX::Zero(0);
    }

    for (auto meas_type : meas_types) {
      switch (meas_type) {
        case LinkMeasurementType::Range:
          rho_ow = GetOneWayRangeMeasurement(epoch_rx, rv_tx, rv_rx, clk_tx, clk_rx, H_ow_range,
                                             hardware_delay, with_noise, with_jacobian);
          z(idx) = rho_ow;
          H_ow_rx.row(idx) = H_ow_range;
          idx++;

        case LinkMeasurementType::RangeRate:
          rho_ow_rate = GetOneWayRangeRateMeasurement(epoch_rx, rv_tx, rv_rx, clk_tx, clk_rx,
                                                      H_ow_rangerate, hardware_delay, with_noise,
                                                      with_jacobian);
          z(idx) = rho_ow_rate;
          H_ow_rx.row(idx) = H_ow_rangerate;
          idx++;

        default: break;
      }
    }

    return z;
  }

  Real LinkMeasurement::GetOneWayRangeMeasurement(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx,
                                                  Vec2 clk_tx, Vec2 clk_rx, MatXd &H_ow_rx,
                                                  Real hardware_delay, bool with_noise,
                                                  bool with_jacobian) {
    // return

    auto func = [epoch_rx, hardware_delay, this](const Vec6 rv_tx_in, const Vec2 clk_tx_in,
                                                 const Vec6 rv_rx_in, const Vec2 clk_rx_in) {
      Real owr = ComputeOneWayRangeLTR(epoch_rx, rv_tx_in, rv_rx_in, clk_tx_in(0), clk_rx_in(0),
                                       linkparams_.tx_center_body, linkparams_.rx_center_body,
                                       linkparams_.is_bodyfixed_tx, linkparams_.is_bodyfixed_rx,
                                       hardware_delay);
      return owr;
    };

    // Without lighttime delay
    // auto func = [rv_tx, clk_tx](const Vec6 rv_rx_in, const Vec2 clk_rx_in) {
    //   Real owr = ComputePseudorange(
    //       rv_tx.head(3), rv_rx_in.head(3), clk_tx(0), clk_rx_in(0), 0.0);
    //   return owr;
    // };

    // break the computational graph relations before taking the jacobian
    Vec6 rv_tx_tmp = rv_tx.cast<double>();
    Vec2 clk_tx_tmp = clk_tx.cast<double>();
    Vec6 rv_rx_tmp = rv_rx.cast<double>();
    Vec2 clk_rx_tmp = clk_rx.cast<double>();
    Real rho_ow = 0.0;

    if (with_jacobian) {
      VecXd H_ow_vec = gradient(func, wrt(rv_tx_tmp, clk_tx_tmp, rv_rx_tmp, clk_rx_tmp),
                                at(rv_tx_tmp, clk_tx_tmp, rv_rx_tmp, clk_rx_tmp), rho_ow);
      // Convert to (1, 8) matrix
      H_ow_rx.row(0) = H_ow_vec.transpose();
    } else {
      rho_ow = func(rv_tx_tmp, clk_tx_tmp, rv_rx_tmp, clk_rx_tmp);
    }

    if (with_noise) {
      double sigma_ow = 0.0;
      // For one-way link, double the range error of the two way link
      if (linkparams_.is_groundstation_rx) {
        sigma_ow = 2
                   * ComputePnRangeErrorCTL(linkparams_.CN0_linear, linkparams_.B_L_chip,
                                            linkparams_.Tc, linkparams_.modulation_type);
      } else {
        sigma_ow = 2
                   * ComputePnRangeErrorOL(linkparams_.CN0_linear, linkparams_.T_I_range,
                                           linkparams_.Tc, linkparams_.modulation_type);
      }

      rho_ow += SampleRandNormal(0.0, sigma_ow, seed_);

      return rho_ow;
    }
  }

  Real LinkMeasurement::GetOneWayRangeRateMeasurement(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx,
                                                      Vec2 clk_tx, Vec2 clk_rx, MatXd &H_ow_rx,
                                                      Real hardware_delay, bool with_noise,
                                                      bool with_jacobian) {
    auto func = [epoch_rx, rv_tx, clk_tx, hardware_delay, this](
                    const Vec6 rv_tx_in, const Vec2 clk_tx_in, const Vec6 rv_rx_in,
                    const Vec2 clk_rx_in) {
      Real owrr = ComputeOneWayRangeRateLTR(
          epoch_rx, rv_tx_in, rv_rx_in, clk_tx(1), clk_rx_in(1), linkparams_.tx_center_body,
          linkparams_.rx_center_body, linkparams_.is_bodyfixed_tx, linkparams_.is_bodyfixed_rx,
          hardware_delay, linkparams_.T_I_doppler);
      return owrr;
    };

    // break the computational graph relations before taking the jacobian
    Vec6 rv_tx_tmp = rv_rx.cast<double>();
    Vec2 clk_tx_tmp = clk_rx.cast<double>();
    Vec6 rv_rx_tmp = rv_rx.cast<double>();
    Vec2 clk_rx_tmp = clk_rx.cast<double>();
    Real rho_ow_rate = 0.0;

    if (with_jacobian) {
      VecXd H_ow_rate_vec = gradient(func, wrt(rv_tx_tmp, clk_tx_tmp, rv_rx_tmp, clk_rx_tmp),
                                     at(rv_tx_tmp, clk_tx_tmp, rv_rx_tmp, clk_rx_tmp), rho_ow_rate);

      // Convert to (1, 8) matrix
      H_ow_rx.row(0) = H_ow_rate_vec.transpose();

    } else {
      rho_ow_rate = func(rv_tx_tmp, clk_tx_tmp, rv_rx_tmp, clk_rx_tmp);
    }

    if (with_noise) {
      double sigma_ow_rate = ComputeRangeRateErrorOneWay(
          linkparams_.B_L_carrier, linkparams_.freq, linkparams_.Tc, linkparams_.T_I_doppler,
          trans_ow_.CN0_linear, linkparams_.sigma_y_1s, linkparams_.modulation_type,
          linkparams_.m_R);

      rho_ow_rate += SampleRandNormal(0.0, sigma_ow_rate, seed_);
    }

    return rho_ow_rate;
  }

  /********************** Two way Link ***************************/

  void LinkMeasurement::GenerateTwoWayLink(Real epoch, std::shared_ptr<Transponder> &tr_receiver,
                                           std::shared_ptr<Transponder> &tr_target,
                                           std::string txrx) {
    SpaceChannel sc = SpaceChannel();
    sc.SetOccultationBodies(occult_bodies_, occult_alt_);

    Real t_tx_u, t_rx_u, t_tx_d, t_rx_d;
    ITransmission trans_d, trans_u;

    // tr_target -> tr_receiver (downlink)
    std::shared_ptr<Transmitter> tx_d = tr_target->GetTransmitter();
    std::shared_ptr<Receiver> rx_d = tr_receiver->GetReceiver();

    // tr_receiver -> tr_target (uplink)
    std::shared_ptr<Transmitter> tx_u = tr_receiver->GetTransmitter();
    std::shared_ptr<Receiver> rx_u = tr_target->GetReceiver();

    Real delay_tx_receiver = hardware_delay_;
    Real delay_rx_receiver = hardware_delay_;
    Real delay_tx_target = hardware_delay_;
    Real delay_rx_target = hardware_delay_;

    if (txrx == "rx") {
      t_rx_d = epoch;
      trans_d = sc.ComputeLinkBudget(tx_d, rx_d, t_rx_d - delay_rx_receiver, "rx");
      t_tx_d = trans_d.t_tx - delay_tx_target;
      t_rx_u = t_tx_d;
      trans_u = sc.ComputeLinkBudget(tx_u, rx_u, t_rx_u - delay_rx_target, "rx");
      t_tx_u = trans_u.t_tx - delay_tx_receiver;

    } else if (txrx == "tx") {
      t_tx_u = epoch;
      trans_u = sc.ComputeLinkBudget(tx_u, rx_u, t_tx_u + delay_tx_receiver, "tx");
      t_rx_u = trans_d.t_rx + delay_rx_target;
      t_tx_d = t_rx_u;
      trans_d = sc.ComputeLinkBudget(tx_d, rx_d, t_tx_d + delay_tx_target, "tx");
      t_rx_d = trans_d.t_rx + delay_rx_receiver;
    }

    // Set time
    epoch_tx_true_ = t_tx_u;
    epoch_rx_true_ = t_rx_d;

    // Set link
    trans_tw_.push_back(trans_u);
    trans_tw_.push_back(trans_d);

    // Register the agents
    agents_.push_back(tr_receiver->GetAgent());
    agents_.push_back(tr_target->GetAgent());

    // Link Parameters -------------------------------------------------
    linkparams_.turnaround_ratio = tr_target->turnaround_ratio;

    SetLinkParams();

    two_way_generated_ = true;
  }

  VecX LinkMeasurement::GetTrueTwoWayLinkMeasurement(std::vector<LinkMeasurementType> meas_types) {
    if (!two_way_generated_) {
      std::cerr << "Error: Two way link not generated" << std::endl;
      return VecX::Zero(0);
    }

    VecX rv_receiver = agents_[0]->GetRvStateAtEpoch(epoch_rx_true_);
    VecX rv_target = agents_[0]->GetRvStateAtEpoch(epoch_rx_true_);

    MatXd H_tw_rx(2, state_size_tw_);  // temporary, won't be used

    return GetTwoWayLinkMeasurement(epoch_rx_true_, rv_receiver, rv_target, H_tw_rx,
                                    hardware_delay_, meas_types, true, false);
  }

  VecX LinkMeasurement::GetTwoWayLinkMeasurement(Real epoch_rx, Vec6 rv_receiver, Vec6 rv_target,
                                                 MatXd H_tw_rx, Real hardware_delay,
                                                 std::vector<LinkMeasurementType> meas_types,
                                                 bool with_noise, bool with_jacobian) {
    Real rho_tw, rho_tw_rate;
    double sigma_tw = 0.0;
    double sigma_tw_rate = 0.0;

    int meas_size = meas_types.size();
    int state_size = state_size_tw_;

    VecX z(meas_size);
    H_tw_rx.resize(meas_size, state_size);

    MatXd H_tw_range = MatXd::Zero(1, state_size);
    MatXd H_tw_rangerate = MatXd::Zero(1, state_size);

    int idx = 0;

    if (meas_size == 0) {
      return VecX::Zero(0);
    }

    for (auto meas_type : meas_types) {
      switch (meas_type) {
        case LinkMeasurementType::Range:
          rho_tw = GetTwoWayRangeMeasurement(epoch_rx, rv_receiver, rv_target, H_tw_range,
                                             hardware_delay, with_noise, with_jacobian);
          z(idx) = rho_tw;
          H_tw_rx.row(idx) = H_tw_range;
          idx++;

        case LinkMeasurementType::RangeRate:
          rho_tw_rate
              = GetTwoWayRangeRateMeasurement(epoch_rx, rv_receiver, rv_target, H_tw_rangerate,
                                              hardware_delay, with_noise, with_jacobian);
          z(idx) = rho_tw_rate;
          H_tw_rx.row(idx) = H_tw_rangerate;
          idx++;

        default: break;
      }
    }

    return z;
  }

  Real LinkMeasurement::GetTwoWayRangeMeasurement(Real epoch_rx, Vec6 rv_receiver, Vec6 rv_target,
                                                  MatXd &H_tw_range, Real hardware_delay,
                                                  bool with_noise, bool with_jacobian) {
    // return

    auto func
        = [epoch_rx, hardware_delay, this](const Vec6 rv_target_in, const Vec6 rv_receiver_in) {
            Real twr = ComputeTwoWayRangeLTR(epoch_rx, rv_target_in, rv_receiver_in,
                                             linkparams_.tx_center_body, linkparams_.rx_center_body,
                                             linkparams_.is_bodyfixed_tx,
                                             linkparams_.is_bodyfixed_rx, hardware_delay, 0.0);
            return twr;
          };

    // break the computational graph relations before taking the jacobian
    Vec6 rv_target_tmp = rv_target.cast<double>();
    Vec6 rv_receiver_tmp = rv_receiver.cast<double>();
    Real rho_tw = 0.0;

    if (with_jacobian) {
      VecXd H_tw_vec = gradient(func, wrt(rv_target_tmp, rv_receiver_tmp),
                                at(rv_target_tmp, rv_receiver_tmp), rho_tw);
      // Convert to (1, state) matrix
      H_tw_range.row(0) = H_tw_vec.transpose();
    } else {
      rho_tw = func(rv_target_tmp, rv_receiver_tmp);
    }

    if (with_noise) {
      double sigma_tw = 0.0;
      if (linkparams_.is_groundstation_rx) {
        sigma_tw = ComputePnRangeErrorCTL(linkparams_.CN0_linear, linkparams_.B_L_chip,
                                          linkparams_.Tc, linkparams_.modulation_type);
      } else {
        sigma_tw = ComputePnRangeErrorOL(linkparams_.CN0_linear, linkparams_.T_I_range,
                                         linkparams_.Tc, linkparams_.modulation_type);
      }

      rho_tw += SampleRandNormal(0.0, sigma_tw, seed_);

      return rho_tw;
    }
  }

  Real LinkMeasurement::GetTwoWayRangeRateMeasurement(Real epoch_rx, Vec6 rv_receiver,
                                                      Vec6 rv_target, MatXd &H_tw_rr,
                                                      Real hardware_delay, bool with_noise,
                                                      bool with_jacobian) {
    // Function to compute the two way range rate
    auto func = [epoch_rx, rv_target, hardware_delay, this](const Vec6 rv_target_in,
                                                            const Vec6 rv_receiver_in) {
      Real owrr = ComputeTwoWayRangeRateLTR(
          epoch_rx, rv_target_in, rv_receiver_in, linkparams_.tx_center_body,
          linkparams_.rx_center_body, linkparams_.is_bodyfixed_tx, linkparams_.is_bodyfixed_rx,
          hardware_delay, linkparams_.T_I_doppler);
      return owrr;
    };

    // break the computational graph relations before taking the jacobian
    Vec6 rv_target_tmp = rv_target.cast<double>();
    Vec6 rv_rx_tmp = rv_receiver.cast<double>();
    Real rho_tw_rate = 0.0;

    if (with_jacobian) {
      VecXd H_tw_rate_vec = gradient(func, wrt(rv_target_tmp, rv_rx_tmp),
                                     at(rv_target_tmp, rv_rx_tmp), rho_tw_rate);

      // Convert to (1, 8) matrix
      H_tw_rr.row(0) = H_tw_rate_vec.transpose();

    } else {
      rho_tw_rate = func(rv_target_tmp, rv_rx_tmp);
    }

    if (with_noise) {
      double sigma_tw_rate = ComputeRangeRateErrorTwoWay(
          linkparams_.B_L_carrier, linkparams_.freq, linkparams_.Tc, linkparams_.T_I_doppler,
          linkparams_.CN0_linear, linkparams_.sigma_y_1s, linkparams_.turnaround_ratio,
          linkparams_.modulation_type, linkparams_.m_R);

      rho_tw_rate += SampleRandNormal(0.0, sigma_tw_rate, seed_);
    }

    return rho_tw_rate;
  }

}  // namespace lupnt
