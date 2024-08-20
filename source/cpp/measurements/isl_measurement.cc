#include "lupnt/measurements/isl_measurement.h"

#include "lupnt/measurements/radio_measurement.h"
#include "lupnt/measurements/space_channel.h"

namespace lupnt {

  /********************** One way ISL ***************************/

  void IslMeasurement::GenerateOneWayLink(Real epoch_rx, std::shared_ptr<Transmitter> &tx,
                                          std::shared_ptr<Receiver> &rx) {
    SpaceChannel sc = SpaceChannel();
    sc.SetOccultationBodies(occult_bodies_, occult_alt_);
    trans_ow_ = sc.ComputeLinkBudget(tx, rx, epoch_rx.val(), "rx");
    one_way_generated_ = true;
  }

  VecX IslMeasurement::GetTrueOneWayISLMeasurement(Real epoch_rx, int agent_id_tx, int agent_id_rx,
                                                   MatXd H_ow_rx,
                                                   std::vector<IslMeasurementType> meas_types) {
    VecX rv_tx = agents_[agent_id_tx]->GetRvStateAtEpoch(epoch_rx);
    VecX rv_rx = agents_[agent_id_rx]->GetRvStateAtEpoch(epoch_rx);
    VecX clk_tx = agents_[agent_id_tx]->GetClockStateVecAtEpoch(epoch_rx);
    VecX clk_rx = agents_[agent_id_rx]->GetClockStateVecAtEpoch(epoch_rx);

    return GetOneWayIslMeasurement(epoch_rx, rv_tx, rv_rx, clk_tx, clk_rx, H_ow_rx, hardware_delay_,
                                   meas_types, true);
  }

  VecX IslMeasurement::GetOneWayIslMeasurement(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx,
                                               Vec2 clk_rx, MatXd H_ow_rx, Real hardware_delay,
                                               std::vector<IslMeasurementType> meas_types,
                                               bool with_noise) {
    Real rho_ow, rho_ow_rate;
    double sigma_ow = 0.0;
    double sigma_ow_rate = 0.0;

    int meas_size = meas_types.size();
    int state_size = 8;

    VecX z(meas_size);
    H_ow_rx.resize(meas_size, state_size);

    int idx = 0;
    Body tx_center_body = CreateDefaultBody(trans_ow_.tx->GetAgent()->GetBodyId());
    Body rx_center_body = CreateDefaultBody(trans_ow_.rx->GetAgent()->GetBodyId());
    bool is_bodyfixed_tx = trans_ow_.tx->GetAgent()->IsBodyFixed();
    bool is_bodyfixed_rx = trans_ow_.rx->GetAgent()->IsBodyFixed();

    if (meas_size == 0) {
      return VecX::Zero(0);
    }

    for (auto meas_type : meas_types) {
      switch (meas_type) {
        case IslMeasurementType::OWR:
          rho_ow = GetOneWayRangeMeasurement(epoch_rx, rv_tx, rv_rx, clk_tx, clk_rx, H_ow_rx,
                                             hardware_delay, &tx_center_body, &rx_center_body,
                                             is_bodyfixed_tx, is_bodyfixed_rx, with_noise);
          z(idx) = rho_ow;
          idx++;
      }
    }

    return z;
  }

  Real IslMeasurement::GetOneWayRangeMeasurement(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx,
                                                 Vec2 clk_rx, MatXd &H_ow_rx, Real hardware_delay,
                                                 Body *tx_center_body, Body *rx_center_body,
                                                 bool is_bodyfixed_tx, bool is_bodyfixed_rx,
                                                 bool with_noise) {
    // return

    auto func = [epoch_rx, rv_tx, clk_tx, tx_center_body, rx_center_body, is_bodyfixed_tx,
                 is_bodyfixed_rx, hardware_delay](Vec6 rv_rx_in, Vec2 clk_rx_in) {
      Real owr = RadioMeasurement::ComputeOneWayRangeLTR(
          epoch_rx, rv_tx, rv_rx_in, clk_tx(0), clk_rx_in(0), tx_center_body, rx_center_body,
          is_bodyfixed_tx, is_bodyfixed_rx, hardware_delay);
      return owr;
    };

    // break the computational graph relations before taking the jacobian
    // Vec6 rv_rx_tmp = rv_rx.cast<double>();
    // Vec2 clk_rx_tmp = clk_rx.cast<double>();

    Vec6 rv_rx_tmp = rv_rx;
    Vec2 clk_rx_tmp = clk_rx;

    Real rho_ow = 0.0;

    // H_ow_rx = jacobian(func, wrt(rv_rx_tmp, clk_rx_tmp),
    //                    at(rv_rx_tmp, clk_rx_tmp), rho_ow);

    if (with_noise) {
      double sigma_ow = RadioMeasurement::ComputePnRangeErrorCTL(
          trans_ow_.CN0_linear, ow_params_.B_L_chip, ow_params_.Tc_chip, ow_params_.carrier_type);
      rho_ow += SampleRandNormal(0.0, sigma_ow, seed_);
    }

    return rho_ow;
  }

  /********************** Two way ISL ***************************/

  void IslMeasurement::GenerateTwoWayLink(Real epoch_rx, std::shared_ptr<Transceiver> &tr1,
                                          std::shared_ptr<Transceiver> &tr2,
                                          std::vector<NaifId> occult_bodies, VecXd occult_alt) {
    SpaceChannel sc = SpaceChannel();
    sc.SetOccultationBodies(occult_bodies, occult_alt);

    // tr2 -> tr1 (downlink)
    std::shared_ptr<Transmitter> tx_d = tr2->GetTransmitter();
    std::shared_ptr<Receiver> rx_d = tr1->GetReceiver();
    ITransmission trans_d = sc.ComputeLinkBudget(tx_d, rx_d, epoch_rx.val(), "rx");

    // tr1 -> tr2 (uplink)
    Real epoch_rx_u = trans_d.t_tx;
    std::shared_ptr<Transmitter> tx_u = tr1->GetTransmitter();
    std::shared_ptr<Receiver> rx_u = tr2->GetReceiver();
    ITransmission trans_u = sc.ComputeLinkBudget(tx_u, rx_u, epoch_rx_u, "rx");

    trans_tw_.push_back(trans_u);
    trans_tw_.push_back(trans_d);

    two_way_generated_ = true;
  }

}  // namespace lupnt
