/**
 * @file gnss_constellation.h
 * @author Stanford NAV LAB
 * @brief GNSS constellation with precomputed ephemerides and transmitter
 *        metadata.
 * @version 0.1
 * @date 2025-06-07
 *
 * @copyright Copyright (c) 2025
 *
 * This class is the C++ counterpart of `pylupnt.measurements.gnss_meas.GNSSMeas`
 * (see `python/pylupnt/measurements/gnss_meas.py`). It is intended to support
 * cislunar / sidelobe GNSS navigation scenarios. Receiver-dependent operations
 * such as light-time iteration, visibility, C/N0 thresholding, and channel
 * selection are intentionally performed by the GNSS measurement model because
 * they depend on the receiver state and receive epoch.
 *
 * Unlike `GNSSMeas`, which propagates / loads GNSS ephemerides itself (via
 * SP3/BRDC loaders), `GnssConstellation` consumes precomputed ephemerides
 * (position/velocity history in ECI) -- either set directly with
 * `SetSatelliteStates`, or loaded from an HDF5 file with `LoadEphemeris`
 * (e.g. produced by a Python preprocessing step using `SP3Loader`/`BRDCLoader`).
 */
#pragma once

#include <filesystem>
#include <functional>
#include <map>
#include <vector>

#include "lupnt/agents/constellation.h"
#include "lupnt/agents/gnss_attitude.h"
#include "lupnt/devices/gnss_device.h"
#include "lupnt/measurements/antenna.h"
#include "lupnt/numerics/cheby_fit.h"

namespace lupnt {

  /// @brief GNSS receiver tracking-loop parameters.
  /// Mirrors `pylupnt.measurements.gnss_meas.GNSSReceiverParam`.
  /// Reference: https://www.mdpi.com/1424-8220/16/3/347
  struct GnssReceiverParams {
    Real Bp = 1.0;   // [Hz] Carrier loop noise bandwidth
    Real T = 20e-3;  // [s] Tracking loop integration time
    Real b = 2.0;    // [-] Front-end bandwidth factor (B_fe = b * R_c)
    Real Bn = 0.7;   // [Hz] Code loop noise bandwidth
    Real Bf = 0.2;   // [Hz] Frequency loop noise bandwidth
    Real D = 0.1;    // [chip] Early-to-late correlator spacing

    // Link-budget terms used by GNSSMeasurements::ComputeCN0 (link-budget / CN0 estimation).
    Real L_ad = 0.6;      // [dB] A/D converter loss
    Real L_pol = 1.0;     // [dB] Polarization loss
    Real L_atm = 0.0;     // [dB] Atmospheric loss
    Real T_eff = 167.98;  // [K] Effective noise temperature
  };

  /// @brief A spherical body that can occlude the line of sight between a
  /// GNSS transmitter and a receiver (e.g. the Earth or the Moon).
  struct GnssOccludingBody {
    Real radius_m = 0.0;             // [m] Body radius
    Vec3 position_m = Vec3::Zero();  // [m] Body center (used when position_provider is null)
    /// Per-epoch position callback; if set, overrides position_m in BuildChannels.
    /// Called with the receive epoch in options_.receive_time_scale; must return [m]
    /// in the same frame as the receiver/transmitter states (options_.frame).
    std::function<Vec3(Real)> position_provider;
  };

  /// @brief Constellation of GNSS satellites with precomputed ephemerides.
  ///
  /// Provides satellite state interpolation plus transmitter antenna/power
  /// metadata. Receiver-dependent measurement construction lives in
  /// `GNSSMeasurements`.
  class GnssConstellation : public Constellation {
  public:
    /// @brief Construct an empty constellation (no PRNs/ephemerides set).
    GnssConstellation() = default;

    /// @brief Construct an empty constellation for a given GNSS system,
    /// naming it after `gnss_const` (e.g. "GPS", "GALILEO").
    /// @param gnss_const GNSS system (GPS, Galileo, BDS, QZSS, ...)
    explicit GnssConstellation(GnssConst gnss_const);

    /// @brief Construct from a YAML config node specifying `gnss_const`
    /// (and optionally `name`).
    ///
    /// Called by the simulation/application builder when a `gnss_const`
    /// section appears in a scenario config; satellite ephemerides must
    /// still be populated afterwards via `SetSatelliteStates`,
    /// `LoadEphemeris`, or `SetupSatelliteStatesFromFiles`.
    /// @param config YAML config node with a `gnss_const` key
    GnssConstellation(Config& config);

    // ---- Setup --------------------------------------------------------------

    /// @brief Set precomputed satellite ephemerides directly from
    /// already-loaded position/velocity histories.
    ///
    /// The primary entry point used by `LoadEphemeris` and
    /// `SetupSatelliteStatesFromFiles`: stores the raw `rv_eci` history per
    /// PRN and fits a piecewise-Chebyshev model (`FitStateHistoryChebyshev`)
    /// to each, which `GetSatelliteStateEci`/`ComputeRange`/`ComputeRangeRate`
    /// then interpolate at arbitrary query epochs.
    /// @param prns PRNs of the satellites (size N)
    /// @param t_tai Epochs of the ephemeris samples, in TAI seconds (size M)
    /// @param rv_eci Per-satellite ECI position/velocity history; each entry
    ///        is an [M x 6] matrix of [r; v] rows (size N)
    void SetSatelliteStates(const std::vector<int>& prns, const VecXd& t_tai,
                            const std::vector<MatXd>& rv_eci);

    /// @brief Load precomputed satellite ephemerides from an HDF5 file.
    /// Expected layout (see also `SaveEphemeris`):
    ///   /prns          [N]   integer PRNs
    ///   /t_tai         [M]   epochs in TAI seconds
    ///   /rv_eci/<prn>  [Mx6] ECI position/velocity history for each PRN
    /// This lets ephemerides be precomputed in Python (e.g. via `SP3Loader` /
    /// `BRDCLoader`, see `GNSSMeas.setup_gnss`) and consumed directly in C++,
    /// avoiding re-implementation of SP3/BRDC parsing.
    void LoadEphemeris(const std::filesystem::path& filepath);

    /// @brief Save the currently-set ephemerides to an HDF5 file using the
    /// layout expected by `LoadEphemeris`.
    void SaveEphemeris(const std::filesystem::path& filepath) const;

    /// @brief Build precomputed satellite ephemerides directly from SP3
    /// precise-ephemeris and ANTEX antenna files (using `Sp3Loader` /
    /// `AntexLoader`, see `lupnt/interfaces/sp3_loader.h` /
    /// `lupnt/interfaces/antex_loader.h`), applying the antenna
    /// phase-center-offset (PCO) correction to the SP3 (center-of-mass) ECEF
    /// positions before converting to ECI -- mirroring the workflow in
    /// `projects/Plasmasphere_Delay_Datagen/generate_antenna_pco.ipynb`:
    ///   pos_corrected_ecef = pos_sp3_ecef + Cijk(t_tai, pos_sp3_ecef) @ pco_neu
    /// (`Cijk` = `AntexLoader::ComputeIjkToEcefRotation`). The (uncorrected)
    /// SP3 velocity is retained -- the PCO is constant in the satellite body
    /// frame, so its contribution to the ECEF velocity is negligible for
    /// link-budget / visibility purposes.
    ///
    /// This is an alternative, "from source files", entry point to
    /// `SetSatelliteStates` / `LoadEphemeris` that avoids a separate Python
    /// preprocessing step.
    /// @param sp3_filepaths SP3 precise-ephemeris file(s) (e.g. consecutive days)
    /// @param antex_filepath ANTEX file providing the satellite antenna PCO
    /// @param t_tai Epochs at which to sample the corrected ephemeris [TAI seconds]
    /// @param freq Frequency whose phase center (PCO) should be used for the
    ///        correction (e.g. `GnssFreq::L1`)
    /// @param prns PRNs to set up; if empty, every PRN of `gnss_const_` found
    ///        in the SP3 file(s) (and present in the ANTEX file) is used
    void SetupSatelliteStatesFromFiles(const std::vector<std::filesystem::path>& sp3_filepaths,
                                       const std::filesystem::path& antex_filepath,
                                       const VecXd& t_tai, GnssFreq freq,
                                       const std::vector<int>& prns = {});

    /// @brief Build a `GnssTransmitter` device for every PRN currently in the
    /// constellation, harvesting its antenna gain pattern(s) and transmit
    /// power(s) (reuses `GnssTransmitter::InitGps` / `InitGalileo` / `InitQzss`,
    /// which load patterns via `gps_table.csv` / antenna pattern files).
    void SetupTransmitters();

    /// @brief Mark the given PRNs as "faulted" so that `IsFaultPrn` reports
    /// them as unavailable, e.g. to simulate a satellite outage for
    /// fault-detection / robustness testing.
    void SetFaultPrns(const std::vector<int>& prns) { fault_prns_ = prns; }
    // ---- Queries --------------------------------------------------------------

    /// @brief Get the number of satellites (PRNs) currently configured.
    int GetNumSatellites() const { return static_cast<int>(prns_.size()); }
    /// @brief Get the list of PRNs currently configured.
    const std::vector<int>& GetPrns() const { return prns_; }
    /// @brief Get the GNSS system (GPS, Galileo, BDS, QZSS, ...) of this
    /// constellation.
    GnssConst GetGnssConst() const { return gnss_const_; }
    /// @brief Get the list of frequencies for which transmitter metadata is
    /// available (set by `SetupTransmitters`, taken from PRN 0).
    const std::vector<GnssFreq>& GetFreqList() const { return freq_list_; }

    /// @brief Check whether `prn` has been marked faulted via
    /// `SetFaultPrns`.
    ///
    /// Called by `GNSSMeasurements::ComputeMeasurements`
    /// (`lupnt/measurements/gnss_measurement.cc`) to skip faulted PRNs when
    /// building the set of visible/usable measurement channels at an epoch.
    /// @param prn PRN to check
    /// @return    True if `prn` is in the fault list
    bool IsFaultPrn(int prn) const;

    /// @brief Interpolated ECI state [r; v] of satellite `prn` at `t_tai`
    /// (linear interpolation over the precomputed ephemeris).
    ///
    /// Evaluates the piecewise-Chebyshev model fitted by
    /// `SetSatelliteStates`; called by `ComputeRange`/`ComputeRangeRate` and
    /// by `GNSSMeasurements` to get a transmitter's state at an
    /// (iteratively-refined) transmit epoch.
    /// @param prn   Satellite PRN
    /// @param t_tai Query epoch [TAI seconds]
    /// @return      ECI state `[r; v]` [m, m/s]
    Vec6 GetSatelliteStateEci(int prn, Real t_tai) const;

    /// @brief Geometric range between satellite `prn` and a receiver position
    /// (ECI) at `t_tai` [m].
    Real ComputeRange(int prn, const Vec3& r_rx_eci, Real t_tai) const;

    /// @brief Range rate (line-of-sight projection of relative velocity)
    /// between satellite `prn` and receiver state (ECI) at `t_tai` [m/s].
    Real ComputeRangeRate(int prn, const Vec6& rv_rx_eci, Real t_tai) const;

    /// @brief Ephemeris epochs in TAI seconds.
    const VecXd& GetEphemerisTimesTai() const { return t_tai_; }

    /// @brief ECI state history [M x 6] for `prn`.
    const MatXd& GetSatelliteStateHistoryEci(int prn) const;

    /// @brief Piecewise-Chebyshev ECI state model [r; v] for `prn`.
    const ChebyshevFitModel& GetSatelliteStateChebyshevEci(int prn) const;

    /// @brief Return true if antenna and transmit-power metadata exist for
    /// `prn` and `freq`.
    bool HasTransmitterInfo(int prn, GnssFreq freq) const;

    /// @brief Transmit antenna pattern for `prn` and `freq`.
    const Antenna& GetTransmitterAntenna(int prn, GnssFreq freq) const;

    /// @brief Transmit power [dB-W] for `prn` and `freq`.
    Real GetTransmitPowerDbw(int prn, GnssFreq freq) const;

    /// @brief Compute the GNSS attitude frame (ex, ey, ez) of a satellite at
    /// position `r_sat_eci`, given the Sun position `r_sun_eci` (same frame).
    /// `ez` points to nadir (toward the central body), `ey` is along the
    /// cross-track direction (toward the Sun), and `ex` completes the right
    /// -handed triad. Mirrors the attitude computation in
    /// `GNSSMeas.setup_measurements`. Delegates to `GnssAttitude::Compute`
    /// (see `lupnt/agents/gnss_attitude.h`).
    static void ComputeAttitude(const Vec3& r_sat_eci, const Vec3& r_sun_eci, Vec3& ex, Vec3& ey,
                                Vec3& ez) {
      GnssAttitude::Compute(r_sat_eci, r_sun_eci, ex, ey, ez);
    }

    /// @brief Compute the GNSS attitude frame (ex, ey, ez) of a satellite at
    /// position `r_sat_eci` with velocity `v_sat_eci`, given the Sun position
    /// `r_sun_eci` (all in the same inertial frame), via the documented
    /// nominal yaw-steering law (`GnssYawSteering::NominalYawAngle`; Eq. 1 of
    /// Cheng et al., 2025). Numerically identical to the Sun-direction-based
    /// `ComputeAttitude` overload above (both describe the same Sun-pointing
    /// nominal attitude), but makes the dependency on the implemented
    /// yaw-steering law explicit. Delegates to `GnssAttitude::Compute`
    /// (see `lupnt/agents/gnss_attitude.h`).
    static void ComputeAttitude(const Vec3& r_sat_eci, const Vec3& v_sat_eci, const Vec3& r_sun_eci,
                                Vec3& ex, Vec3& ey, Vec3& ez) {
      GnssAttitude::Compute(r_sat_eci, v_sat_eci, r_sun_eci, ex, ey, ez);
    }

    /// @brief Pseudorange (DLL/code-tracking) measurement noise std [m].
    /// Mirrors `GNSSMeas.compute_gnss_pseudorange_noise`.
    Real ComputeSigmaRange(Real cn0_dbhz, GnssFreq freq) const;

    /// @brief Pseudorange-rate (FLL/frequency-tracking) measurement noise
    /// std [m/s]. Mirrors `GNSSMeas.compute_gnss_pseudorangerate_noise`.
    Real ComputeSigmaRangeRate(Real cn0_dbhz, GnssFreq freq) const;

    /// @brief Carrier-phase (PLL/phase-tracking) measurement noise std [m].
    /// Mirrors `GNSSMeas.compute_gnss_carrier_phase_noise`.
    Real ComputeSigmaCarrierPhase(Real cn0_dbhz, GnssFreq freq) const;

    // ---- Constellation interface ------------------------------------------
    // Note: `Constellation::Setup`/`Step` are not virtual; `GnssConstellation`
    // does not populate the base `satellites_` container (its GNSS satellites
    // are represented by precomputed ephemerides + lightweight `GnssTransmitter`
    // metadata rather than full dynamically-propagated `Agent`s), so these
    // hide (rather than override) the base no-op implementations.
    void Setup();
    void Step(Real t);

  private:
    GnssConst gnss_const_ = GnssConst::GPS;
    std::vector<int> prns_;
    std::vector<int> fault_prns_;

    VecXd t_tai_;                                    // [M] ephemeris epochs (TAI seconds)
    std::map<int, MatXd> rv_eci_;                    // PRN -> [M x 6] ECI state history
    std::map<int, ChebyshevFitModel> rv_eci_cheby_;  // PRN -> fitted ECI state model

    std::map<int, std::map<GnssFreq, Antenna>> antennas_;  // PRN -> freq -> antenna pattern
    std::map<int, std::map<GnssFreq, Real>> P_tx_;         // PRN -> freq -> transmit power [dB-W]
    std::vector<GnssFreq> freq_list_;                      // Frequencies available (from PRN 0)

    GnssReceiverParams rx_params_;
  };

}  // namespace lupnt
