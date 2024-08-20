/**
 * @file intersatellite_link.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-04-04
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include <vector>

#include "lupnt/agents/agent.h"
#include "lupnt/core/constants.h"
#include "radio_measurement.h"
#include "transmission.h"

namespace lupnt {

enum class IslMeasurementType {
  OWR,   // One-way range
  OWRR,  // One-way range rate
  TWR,   // Two-way range
  TWRR,  // Two-way range rate
  DWR,   // Dual one-way range
  DWRR,  // Dual one-way range rate
};

struct LinkParams {
  double freq;               // carrier frequency [Hz]
  CarrierType carrier_type;  // carrier type
  double B_L_chip;           // tracking loop noise bandwidth
  double Tc;                 // chip duration
  double B_L_carrier;        // carrier loop noise bandwidth
  double sigma_y_1s;         // one-sigma range noise [m]
  double m_R;                // modulation index
  double T_I;                // Doppler integration time
  double G;                  //

  Real hardware_delay = 1e-6;  // hardware delay [s]
};

class IslMeasurement {
 private:
  // True Received Epoch
  Real epoch_rx_true_;

  // Transmission (map to satellite index)
  ITransmission trans_ow_;                // One-way transmission
  std::vector<ITransmission> trans_tw_;   // Two-way transmission
  std::vector<ITransmission> trans_dow_;  // Dual one-way transmission

  // Parameters
  LinkParams linkparams_;

  // Flag if the link is generated
  bool one_way_generated_ = false;
  bool two_way_generated_ = false;
  bool dual_one_way_generated_ = false;

  // Agents
  std::vector<std::shared_ptr<Agent>> agents_;

  // Random seed
  int seed_ = 0;

  // Occultation bodies
  std::vector<NaifId> occult_bodies_;
  VecXd occult_alt_;

 public:
  /**
   * @brief Construct a new Isl Measurement object (one way)
   *
   * @param epoch_rx_true  true received epoch
   * @param tx  transmitter
   * @param rx  receiver
   * @param occult_bodies  occulting bodies
   * @param occult_alt  occultation altitude
   * @param link_type  link type
   */
  IslMeasurement(Real epoch_rx_true, std::shared_ptr<Transmitter> &tx,
                 std::shared_ptr<Receiver> &rx,
                 std::vector<NaifId> occult_bodies, VecXd occult_alt,
                 std::string link_type);

  /**
   * @brief Construct a new Isl Measurement object (two way)
   *
   * @param epoch_rx_true
   * @param target
   * @param rx
   * @param occult_bodies
   * @param occult_alt
   * @param link_type
   */
  IslMeasurement(Real epoch_rx_true, std::shared_ptr<Transceiver> &target,
                 std::shared_ptr<Transceiver> &rx,
                 std::vector<NaifId> occult_bodies, VecXd occult_alt,
                 std::string link_type);

  /********************** One way ISL ***************************/

  /**
   * @brief Get the One Way Link object
   *
   * @param epoch_rx        epoch of the receiver
   * @param tx              transmitter
   * @param rx              receiver
   * @param occult_bodies   occulting bodies
   * @param occult_alt      occultation altitude
   * @return ITransmission
   */
  void GenerateOneWayLink(Real epoch_rx, std::shared_ptr<Transmitter> &tx,
                          std::shared_ptr<Receiver> &rx);

  VecX GetTrueOneWayISLMeasurement(Real epoch_rx, MatXd H_ow_rx,
                                   std::vector<IslMeasurementType> meas_types);

  VecX GetOneWayIslMeasurement(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx,
                               Vec2 clk_tx, Vec2 clk_rx, MatXd H_ow_rx,
                               Real hardware_delay,
                               std::vector<IslMeasurementType> meas_types,
                               bool with_noise);

  Real GetOneWayRangeMeasurement(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx,
                                 Vec2 clk_tx, Vec2 clk_rx, MatXd &H_ow_rx,
                                 Real hardware_delay, Body *tx_center_body,
                                 Body *rx_center_body, bool is_bodyfixed_tx,
                                 bool is_bodyfixed_rx, bool with_noise);

  Real GetOneWayRangeRateMeasurement(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx,
                                     Vec2 clk_tx, Vec2 clk_rx, MatXd &H_ow_rx,
                                     Real hardware_delay, Body *tx_center_body,
                                     Body *rx_center_body, bool is_bodyfixed_tx,
                                     bool is_bodyfixed_rx, bool with_noise);

  /********************** Two way ISL ***************************/

  /**
   * @brief Generate Two way link (receiver->target->receiver)
   *
   * @param epoch_rx       epoch of the receiver
   * @param tr_receiver    transceiver 1 (receiver)
   * @param tr_target      transceiver 2 (target)
   * @param occult_bodies  occulting bodies
   * @param occult_alt     occultation altitude
   * @return std::vector<ITransmission>
   */
  void GenerateTwoWayLink(Real epoch_rx,
                          std::shared_ptr<Transceiver> &tr_receiver,
                          std::shared_ptr<Transceiver> &tr_target);
};

}  // namespace lupnt
