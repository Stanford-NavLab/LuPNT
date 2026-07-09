#pragma once

#include "lupnt/devices/comm_devices.h"
#include "lupnt/measurements/antenna.h"

namespace lupnt {

  enum class GnssConst {
    GPS,      // Global Positioning System
    GLONASS,  // Global Navigation Satellite System
    GALILEO,  // European Global Navigation Satellite System
    BEIDOU,   // Chinese Global Navigation Satellite System
    QZSS,     // Quasi-Zenith Satellite System
  };

  enum class GnssFreq { L1, L2, L5, E1, E6, E5, E5a, E5b };

  extern const std::map<GnssFreq, Real> GNSS_FREQ_MAP;
  extern const std::map<GnssFreq, Real> GNSS_RC_MAP;

  /// @brief GNSS navigation-signal transmitter device: holds per-frequency
  /// transmit power, chip rate, and antenna gain pattern for a specific GNSS
  /// satellite (constellation + PRN).
  ///
  /// Attached to a `GnssConstellation` agent (or constructed standalone, as
  /// in `GnssConstellation::SetupTransmitters`) to provide the
  /// `GetTransmitPower`/`GetAntennas`/`GetFreqList` data consumed by the GNSS
  /// link-budget and pseudorange/carrier-phase measurement models in
  /// `measurements/` to compute received signal power and antenna gain
  /// corrections.
  class GnssTransmitter : public Transmitter {
  public:
    /// @brief Construct directly from a GNSS constellation type and PRN,
    /// dispatching to the matching `Init*` (`InitGps`/`InitGlonass`/
    /// `InitGalileo`/`InitBeidou`/`InitQzss`) to load transmit power(s) and
    /// antenna pattern(s) for that satellite.
    GnssTransmitter(GnssConst gnss_const, int prn);

    /// @brief Construct from a YAML config node specifying `gnss_const` and
    /// `prn`, then dispatch to the matching `Init*` as in the
    /// `(GnssConst, int)` constructor.
    GnssTransmitter(Config& config);

    /// @brief Get the per-frequency transmit power [dB-W] for this satellite.
    std::map<GnssFreq, Real> GetTransmitPower() { return P_tx_; }
    /// @brief Get the per-frequency code chip rate [Hz] for this satellite
    /// (currently unpopulated by any `Init*`; reserved for future use).
    std::map<GnssFreq, Real> GetChipRate() { return R_c_; }
    /// @brief Get the per-frequency antenna gain-pattern model (`Antenna`)
    /// for this satellite.
    std::map<GnssFreq, Antenna> GetAntennas() { return antennas_; }
    /// @brief Get the list of GNSS signal frequencies broadcast by this
    /// satellite (e.g. `{L1, L2, L5}`).
    std::vector<GnssFreq> GetFreqList() { return freq_list_; }

  protected:
    GnssConst gnss_const_;                  // [-] Type of the GNSS system
    int prn_;                               // [-] PRN of the transmitter satellite
    std::vector<GnssFreq> freq_list_;       // [-] List of frequencies (by signal names)
    std::map<GnssFreq, Antenna> antennas_;  // [-] List of antennas
    std::map<GnssFreq, Real> P_tx_;         // [dB-W] Transmitter power
    std::map<GnssFreq, Real> R_c_;          // [m] Chip rate

    /// @brief Initialize `freq_list_`, `P_tx_`, and `antennas_` for a GPS
    /// satellite by looking up `prn_` in the bundled `gps_table.csv` (block
    /// type -> frequencies, transmit power, and antenna pattern name).
    /// Aborts (`LUPNT_CHECK`) if `prn_` is not found in the table.
    void InitGps();

    /// @brief GLONASS initialization; not yet implemented (always aborts via
    /// `LUPNT_CHECK`).
    void InitGlonass();

    /// @brief Initialize `freq_list_`, `antennas_`, and `P_tx_` for a Galileo
    /// satellite (frequencies E1/E5a/E5b/E6, fixed transmit power, and
    /// per-frequency antenna patterns named `"Galileo_<freq>.txt"`).
    void InitGalileo();

    /// @brief BeiDou initialization; not yet implemented (always aborts via
    /// `LUPNT_CHECK`).
    void InitBeidou();

    /// @brief Initialize `freq_list_`, `antennas_`, and `P_tx_` for a QZSS
    /// satellite based on `prn_` (PRNs 1-4 broadcast L1/L2/L5, others
    /// L1/L5), with per-frequency antenna patterns named
    /// `"QZSS_<id>_<freq>.txt"`.
    void InitQzss();
  };

  /// @brief GNSS receiver device.
  ///
  /// Thin `Receiver` specialization (currently adds no behavior beyond the
  /// base class) attached to a satellite/rover/ground-station `Agent` acting
  /// as a GNSS user receiver in `measurements/` GNSS measurement models.
  class GnssReceiver : public Receiver {
  public:
    GnssReceiver() = default;
    /// @brief Construct a GNSS receiver from a YAML config node (see
    /// `Receiver::Receiver(Config&)`).
    GnssReceiver(Config& config);
  };

}  // namespace lupnt
