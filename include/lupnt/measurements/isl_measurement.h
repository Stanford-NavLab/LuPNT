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

  struct OneWayLinkParams {
    CarrierType carrier_type;  // carrier type
    double B_L_chip;           // tracking loop noise bandwidth
    double Tc_chip;            // chip duration
    double m_R;                // modulation index
  };

  struct TwoWayLinkParams {
    CarrierType carrier_type;
    double B_L_chip;
    double Tc_chip;
    double m_R;
  };

  struct DualOneWayLinkParams {
    CarrierType carrier_type;
    double B_L_chip;
    double Tc_chip;
    double m_R;
  };

  class IslMeasurement {
  private:
    // Register set of agents
    std::vector<std::shared_ptr<Agent>> agents_;  // set of agents

    // Transmission (map to satellite index)
    ITransmission trans_ow_;                // One-way transmission
    std::vector<ITransmission> trans_tw_;   // Two-way transmission
    std::vector<ITransmission> trans_dow_;  // Dual one-way transmission

    std::pair<int, int> agent_pair_ow_;
    std::pair<int, int> agent_pair_tw_;
    std::pair<int, int> agent_pair_dow_;

    // Parameters
    OneWayLinkParams ow_params_;
    TwoWayLinkParams tw_params_;
    DualOneWayLinkParams dow_params_;

    // Flag if the link is generated
    bool one_way_generated_ = false;
    bool two_way_generated_ = false;
    bool dual_one_way_generated_ = false;

    // Random seed
    int seed_ = 0;

    // delay in two-way ranging [s]
    double hardware_delay_ = 1e-6;

    // Occultation bodies
    std::vector<NaifId> occult_bodies_;
    VecXd occult_alt_;

  public:
    IslMeasurement(std::vector<NaifId> occult_bodies, VecXd occult_alt) {
      occult_bodies_ = occult_bodies;
      occult_alt_ = occult_alt;
    }

    void AddAgent(std::shared_ptr<Agent> agent) { agents_.push_back(agent); }
    void SetHardwareDelay(double hardware_delay) { hardware_delay_ = hardware_delay; }

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

    VecX GetTrueOneWayISLMeasurement(Real epoch_rx, int agent_id_tx, int agent_id_rx, MatXd H_ow_rx,
                                     std::vector<IslMeasurementType> meas_types);

    VecX GetOneWayIslMeasurement(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx, Vec2 clk_rx,
                                 MatXd H_ow_rx, Real hardware_delay,
                                 std::vector<IslMeasurementType> meas_types, bool with_noise);

    Real GetOneWayRangeMeasurement(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx, Vec2 clk_rx,
                                   MatXd &H_ow_rx, Real hardware_delay, Body *tx_center_body,
                                   Body *rx_center_body, bool is_bodyfixed_tx, bool is_bodyfixed_rx,
                                   bool with_noise);

    /********************** Two way ISL ***************************/

    /**
     * @brief Generate Two way link
     *
     * @param epoch_rx       epoch of the receiver
     * @param tr1            transceiver 1
     * @param tr2            transceiver 2
     * @param occult_bodies  occulting bodies
     * @param occult_alt     occultation altitude
     * @return std::vector<ITransmission>
     */
    void GenerateTwoWayLink(Real epoch_rx, std::shared_ptr<Transceiver> &tr1,
                            std::shared_ptr<Transceiver> &tr2, std::vector<NaifId> occult_bodies,
                            VecXd occult_alt);
  };

}  // namespace lupnt
