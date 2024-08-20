#include "lupnt/measurements/isl_measurement.h"

#include "lupnt/core/constants.h"
#include "lupnt/measurements/radio_measurement.h"
#include "lupnt/measurements/space_channel.h"

namespace lupnt {

IslMeasurement::IslMeasurement(Real epoch_rx_true,
                               std::shared_ptr<Transmitter> &tx,
                               std::shared_ptr<Receiver> &rx,
                               std::vector<NaifId> occult_bodies,
                               VecXd occult_alt,
                               std::string link_type = "two_way")
    : epoch_rx_true_(epoch_rx_true),
      occult_bodies_(occult_bodies),
      occult_alt_(occult_alt) {
  if (link_type == "one_way") {
    GenerateOneWayLink(epoch_rx_true, tx, rx);
  } else if (link_type == "two_way") {
    std::cerr << "Error: Invalid link type, pass transceivers for two-way link"
              << std::endl;
  }

  // Register the agents
  agents_.push_back(tx->GetAgent());
  agents_.push_back(rx->GetAgent());
}

IslMeasurement::IslMeasurement(Real epoch_rx_true,
                               std::shared_ptr<Transceiver> &target,
                               std::shared_ptr<Transceiver> &receiver,
                               std::vector<NaifId> occult_bodies,
                               VecXd occult_alt, std::string link_type)
    : epoch_rx_true_(epoch_rx_true),
      occult_bodies_(occult_bodies),
      occult_alt_(occult_alt) {
  if (link_type == "two_way") {
    GenerateTwoWayLink(epoch_rx_true, receiver, target);
  } else if (link_type == "one_way") {
    std::cerr << "Error: Invalid link type, pass transmitter and receiver for "
                 "one-way link"
              << std::endl;
  }

  // Register the agents
  agents_.push_back(target->GetAgent());
  agents_.push_back(receiver->GetAgent());
}

/********************** One way ISL ***************************/

void IslMeasurement::GenerateOneWayLink(Real epoch_rx,
                                        std::shared_ptr<Transmitter> &tx,
                                        std::shared_ptr<Receiver> &rx) {
  SpaceChannel sc = SpaceChannel();
  sc.SetOccultationBodies(occult_bodies_, occult_alt_);
  trans_ow_ = sc.ComputeLinkBudget(tx, rx, epoch_rx.val(), "rx");
  one_way_generated_ = true;
}

VecX IslMeasurement::GetTrueOneWayISLMeasurement(
    Real epoch_rx, MatXd H_ow_rx, std::vector<IslMeasurementType> meas_types) {
  if (!one_way_generated_) {
    std::cerr << "Error: One way link not generated" << std::endl;
    return VecX::Zero(0);
  }

  VecX rv_tx = agents_[0]->GetRvStateAtEpoch(epoch_rx);
  VecX rv_rx = agents_[0]->GetRvStateAtEpoch(epoch_rx);
  VecX clk_tx = agents_[1]->GetClockStateVecAtEpoch(epoch_rx);
  VecX clk_rx = agents_[1]->GetClockStateVecAtEpoch(epoch_rx);

  return GetOneWayIslMeasurement(epoch_rx, rv_tx, rv_rx, clk_tx, clk_rx,
                                 H_ow_rx, linkparams_.hardware_delay,
                                 meas_types, true);
}

VecX IslMeasurement::GetOneWayIslMeasurement(
    Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx, Vec2 clk_rx,
    MatXd H_ow_rx, Real hardware_delay,
    std::vector<IslMeasurementType> meas_types, bool with_noise) {
  Real rho_ow, rho_ow_rate;
  double sigma_ow = 0.0;
  double sigma_ow_rate = 0.0;

  int meas_size = meas_types.size();
  int state_size = 8;

  VecX z(meas_size);
  H_ow_rx.resize(meas_size, state_size);

  int idx = 0;
  Body tx_center_body =
      CreateDefaultBody(trans_ow_.tx->GetAgent()->GetBodyId());
  Body rx_center_body =
      CreateDefaultBody(trans_ow_.rx->GetAgent()->GetBodyId());
  bool is_bodyfixed_tx = trans_ow_.tx->GetAgent()->IsBodyFixed();
  bool is_bodyfixed_rx = trans_ow_.rx->GetAgent()->IsBodyFixed();

  if (meas_size == 0) {
    return VecX::Zero(0);
  }

  for (auto meas_type : meas_types) {
    switch (meas_type) {
      case IslMeasurementType::OWR:
        rho_ow = GetOneWayRangeMeasurement(
            epoch_rx, rv_tx, rv_rx, clk_tx, clk_rx, H_ow_rx, hardware_delay,
            &tx_center_body, &rx_center_body, is_bodyfixed_tx, is_bodyfixed_rx,
            with_noise);
        z(idx) = rho_ow;
        idx++;

      case IslMeasurementType::OWRR:
        rho_ow_rate = GetOneWayRangeRateMeasurement(
            epoch_rx, rv_tx, rv_rx, clk_tx, clk_rx, H_ow_rx, hardware_delay,
            &tx_center_body, &rx_center_body, is_bodyfixed_tx, is_bodyfixed_rx,
            with_noise);
        z(idx) = rho_ow_rate;
        idx++;
    }
  }

  return z;
}

Real IslMeasurement::GetOneWayRangeMeasurement(
    Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx, Vec2 clk_rx,
    MatXd &H_ow_rx, Real hardware_delay, Body *tx_center_body,
    Body *rx_center_body, bool is_bodyfixed_tx, bool is_bodyfixed_rx,
    bool with_noise) {
  // return

  auto func = [epoch_rx, rv_tx, clk_tx, tx_center_body, rx_center_body,
               is_bodyfixed_tx, is_bodyfixed_rx,
               hardware_delay](const Vec6 rv_rx_in, const Vec2 clk_rx_in) {
    Real owr = RadioMeasurement::ComputeOneWayRangeLTR(
        epoch_rx, rv_tx, rv_rx_in, clk_tx(0), clk_rx_in(0), tx_center_body,
        rx_center_body, is_bodyfixed_tx, is_bodyfixed_rx, hardware_delay);
    return owr;
  };

  // Without lighttime delay
  // auto func = [rv_tx, clk_tx](const Vec6 rv_rx_in, const Vec2 clk_rx_in) {
  //   Real owr = RadioMeasurement::ComputePseudorange(
  //       rv_tx.head(3), rv_rx_in.head(3), clk_tx(0), clk_rx_in(0), 0.0);
  //   return owr;
  // };

  // break the computational graph relations before taking the jacobian
  Vec6 rv_rx_tmp = rv_rx.cast<double>();
  Vec2 clk_rx_tmp = clk_rx.cast<double>();
  Real rho_ow = 0.0;
  VecXd H_ow_vec = gradient(func, wrt(rv_rx_tmp, clk_rx_tmp),
                            at(rv_rx_tmp, clk_rx_tmp), rho_ow);

  // Convert to (1, 8) matrix
  H_ow_rx.row(0) = H_ow_vec.transpose();

  if (with_noise) {
    double sigma_ow = RadioMeasurement::ComputePnRangeErrorCTL(
        trans_ow_.CN0_linear, linkparams_.B_L_chip, linkparams_.Tc,
        linkparams_.carrier_type);
    rho_ow += SampleRandNormal(0.0, sigma_ow, seed_);
  }

  return rho_ow;
}

Real IslMeasurement::GetOneWayRangeRateMeasurement(
    Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx, Vec2 clk_rx,
    MatXd &H_ow_rx, Real hardware_delay, Body *tx_center_body,
    Body *rx_center_body, bool is_bodyfixed_tx, bool is_bodyfixed_rx,
    bool with_noise) {
  double T_I = linkparams_.T_I;

  auto func = [epoch_rx, rv_tx, clk_tx, tx_center_body, rx_center_body,
               is_bodyfixed_tx, is_bodyfixed_rx, T_I,
               hardware_delay](const Vec6 rv_rx_in, const Vec2 clk_rx_in) {
    Real owrr = RadioMeasurement::ComputeOneWayRangeRateLTR(
        epoch_rx, rv_tx, rv_rx_in, clk_tx(1), clk_rx_in(1), tx_center_body,
        rx_center_body, is_bodyfixed_tx, is_bodyfixed_rx, hardware_delay, T_I);
    return owrr;
  };

  // break the computational graph relations before taking the jacobian
  Vec6 rv_rx_tmp = rv_rx.cast<double>();
  Vec2 clk_rx_tmp = clk_rx.cast<double>();
  Real rho_ow_rate = 0.0;
  VecXd H_ow_rate_vec = gradient(func, wrt(rv_rx_tmp, clk_rx_tmp),
                                 at(rv_rx_tmp, clk_rx_tmp), rho_ow_rate);

  // Convert to (1, 8) matrix
  H_ow_rx.row(1) = H_ow_rate_vec.transpose();

  if (with_noise) {
    double sigma_ow_rate = RadioMeasurement::ComputeRangeRateErrorOneWay(
        linkparams_.B_L_carrier, linkparams_.freq, linkparams_.Tc,
        linkparams_.T_I, trans_ow_.CN0_linear, linkparams_.sigma_y_1s,
        linkparams_.carrier_type, linkparams_.m_R);

    rho_ow_rate += SampleRandNormal(0.0, sigma_ow_rate, seed_);
  }

  return rho_ow_rate;
}

/********************** Two way ISL ***************************/

void IslMeasurement::GenerateTwoWayLink(
    Real epoch_rx, std::shared_ptr<Transceiver> &tr_receiver,
    std::shared_ptr<Transceiver> &tr_target) {
  SpaceChannel sc = SpaceChannel();
  sc.SetOccultationBodies(occult_bodies_, occult_alt_);

  // tr_target -> tr_receiver (downlink)
  std::shared_ptr<Transmitter> tx_d = tr_target->GetTransmitter();
  std::shared_ptr<Receiver> rx_d = tr_receiver->GetReceiver();
  ITransmission trans_d =
      sc.ComputeLinkBudget(tx_d, rx_d, epoch_rx.val(), "rx");

  // tr_receiver -> tr_target (uplink)
  Real epoch_rx_u = trans_d.t_tx;
  std::shared_ptr<Transmitter> tx_u = tr_receiver->GetTransmitter();
  std::shared_ptr<Receiver> rx_u = tr_target->GetReceiver();
  ITransmission trans_u = sc.ComputeLinkBudget(tx_u, rx_u, epoch_rx_u, "rx");

  trans_tw_.push_back(trans_u);
  trans_tw_.push_back(trans_d);

  two_way_generated_ = true;
}

}  // namespace lupnt
