/**
 * @file link_measurement.h
 * @author Stanford NAV LAB
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

  enum class LinkMeasurementType { Range, RangeRate, Bearing };

  struct LinkParams {
    // Receiver Parameters
    double freq;                 // carrier frequency [Hz]
    Modulation modulation_type;  // carrier type
    double B_L_chip;             // tracking loop noise bandwidth
    double Tc;                   // chip duration
    double B_L_carrier;          // carrier loop noise bandwidth
    double sigma_y_1s;           // one-sigma range noise [m]
    double m_R;                  // modulation index
    double T_I_doppler;          // Doppler integration time
    double T_I_range;            // range integration time (for open loop)
    double turnaround_ratio;     // Transponder turnaround ratio

    // Agent Parameters
    BodyData tx_center_body;   // center body of the transmitter (target)
    BodyData rx_center_body;   // center body of the receiver
    bool is_bodyfixed_tx;      // is the transmitter (target) body fixed
    bool is_bodyfixed_rx;      // is the receiver body fixed
    bool is_groundstation_rx;  // is the receiver a ground station
    bool is_groundstation_tx;  // is the transmitter (target) a ground station

    double CN0_linear;  // Carrier-to-noise density [dB-Hz]
  };

  class LinkMeasurement {
  private:
    // True epochs
    Real epoch_tx_true_ = 0.0;  // True transmitted epoch (The acutal signal is
                                // tranmistted at epoch_tx_true + hardware delay)
    Real epoch_rx_true_ = 0.0;  // True received epoch (The acutal signal is
                                // received at epoch_rx_true - hardware delay)

    // Transmission
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

    // state size
    int state_size_ow_ = 16;  // One way link state size (target rv + clock,  receiver rv + clock)
    int state_size_tw_ = 12;  // Two way link state size (target rv + receiver rv)

    // Occultation bodies
    std::vector<NaifId> occult_bodies_;
    VecXd occult_alt_;

    // Use Fixed error for the measurements
    bool use_fixed_error_ = false;
    double range_sigma_fixed_ = 0.0;
    double range_rate_sigma_fixed_ = 0.0;

    // Hardware delay (same for transmitter and receiver)
    Real hardware_delay_ = 1e-9;  // hardware delay [s]

  public:
    /**
     * @brief Construct a new Isl Measurement object using transmitters and
     * receivers
     *
     * @param occult_bodies  occulting bodies
     * @param occult_alt  occultation altitude
     * @param hardware_delay  hardware delay
     * @param link_type  link type  (e.g. "one-way", "two-way", "dual-one-way")
     */
    LinkMeasurement(std::vector<NaifId> occult_bodies, VecXd occult_alt, Real hardware_delay);

    /********************** Utils  *********************************/
    void SetLinkParams();
    inline void SetSeed(int seed) { seed_ = seed; }
    inline void SetFixedRangeError(double range_sigma) { range_sigma_fixed_ = range_sigma; }
    inline void SetFixedRangeRateError(double range_rate_sigma) {
      range_rate_sigma_fixed_ = range_rate_sigma;
    }
    inline void UseFixedError() { use_fixed_error_ = true; }
    inline void DisableFixedError() { use_fixed_error_ = false; }

    inline int GetSeed() const { return seed_; }
    inline LinkParams GetLinkParams() const { return linkparams_; }
    inline Real GetTxEpoch() const { return epoch_tx_true_; }
    inline Real GetRxEpoch() const { return epoch_rx_true_; }

    void Reset() {
      one_way_generated_ = false;
      two_way_generated_ = false;
      dual_one_way_generated_ = false;
      // reset transmission (one-way)
      trans_ow_ = ITransmission();
      trans_tw_.clear();
      trans_dow_.clear();
    }

    /********************** One way Link ***************************/

    /**
     * @brief Get the One Way Link object
     *
     * @param epoch_rx        epoch of the receiver or transmitter
     * @param tx              transmitter
     * @param rx              receiver
     * @param fixed_txrx      fixed transmitter or receiver (tx or rx)
     * @return ITransmission
     */
    void GenerateOneWayLink(Real epoch, std::shared_ptr<Transmitter> &tx,
                            std::shared_ptr<Receiver> &rx, std::string fixed_txrx);

    void GenerateOneWayLinkAtRxEpoch(Real epoch, std::shared_ptr<Transmitter> &tx,
                                     std::shared_ptr<Receiver> &rx) {
      GenerateOneWayLink(epoch, tx, rx, "rx");
    };

    void GenerateOneWayLinkAtTxEpoch(Real epoch, std::shared_ptr<Transmitter> &tx,
                                     std::shared_ptr<Receiver> &rx) {
      GenerateOneWayLink(epoch, tx, rx, "tx");
    };

    VecX GetTrueOneWayLinkMeasurement(std::vector<LinkMeasurementType> meas_types);

    VecX GetOneWayLinkMeasurement(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx, Vec2 clk_rx,
                                  MatXd H_ow_rx, Real hardware_delay,
                                  std::vector<LinkMeasurementType> meas_types, bool with_noise,
                                  bool with_jacobian);

    Real GetOneWayRangeMeasurement(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx, Vec2 clk_rx,
                                   MatXd &H_ow_rx, Real hardware_delay, bool with_noise,
                                   bool with_jacobian);

    Real GetOneWayRangeRateMeasurement(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx, Vec2 clk_tx,
                                       Vec2 clk_rx, MatXd &H_ow_rx, Real hardware_delay,
                                       bool with_noise, bool with_jacobian);

    /********************** Two way Link ***************************/

    /**
     * @brief Generate Two way link (receiver->target->receiver)
     *
     * @param epoch       epoch of the receiver or transmitter
     * @param tr_receiver    Transponder (receiver)
     * @param tr_target      Transponder (target)
     * @param txrx_fixed     fixed time of transmitter or receiver (tx or rx)
     * @return std::vector<ITransmission>
     */
    void GenerateTwoWayLink(Real epoch, std::shared_ptr<Transponder> &tr_receiver,
                            std::shared_ptr<Transponder> &tr_target, std::string txrx_fixed);

    void GenerateTwoWayLinkAtRxEpoch(Real epoch, std::shared_ptr<Transponder> &tr_receiver,
                                     std::shared_ptr<Transponder> &tr_target) {
      GenerateTwoWayLink(epoch, tr_receiver, tr_target, "rx");
    };

    void GenerateTwoWayLinkAtTxEpoch(Real epoch, std::shared_ptr<Transponder> &tr_receiver,
                                     std::shared_ptr<Transponder> &tr_target) {
      GenerateTwoWayLink(epoch, tr_receiver, tr_target, "tx");
    };

    VecX GetTrueTwoWayLinkMeasurement(std::vector<LinkMeasurementType> meas_types);

    VecX GetTwoWayLinkMeasurement(Real epoch_rx, Vec6 rv_tx, Vec6 rv_rx, MatXd H_tw_rx,
                                  Real hardware_delay, std::vector<LinkMeasurementType> meas_types,
                                  bool with_noise, bool with_jacobian);

    Real GetTwoWayRangeMeasurement(Real epoch_rx, Vec6 rv_receiver, Vec6 rv_target,
                                   MatXd &H_tw_range, Real hardware_delay, bool with_noise,
                                   bool with_jacobian);

    Real GetTwoWayRangeRateMeasurement(Real epoch_rx, Vec6 rv_receiver, Vec6 rv_target,
                                       MatXd &H_tw_rr, Real hardware_delay, bool with_noise,
                                       bool with_jacobian);
  };

}  // namespace lupnt
