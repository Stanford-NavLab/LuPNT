#pragma once

#include "lupnt/core/definitions.h"

namespace lupnt {

  struct DllParams {
    Real B_dll;  // [Hz] Bandwidth
    Real T_i;    // [s] Coherent integration time
    Real d;      // [chip] Chip duration
    Real Tc;     // [s] Spread code period
    Real B_fe;   // [Hz] Front-end bandwidth
  };

  struct PllParams {
    Real B_pll;  // [Hz] Bandwidth
    Real T_i;    // [s] Coherent integration time
  };

  struct FllParams {
    Real B_fll;      // [Hz] Bandwidth
    Real T_i;        // [s] Coherent integration time
    Real CN0_F_fll;  // [dB-Hz] CN0 limit for low F
  };

  /// @brief Compute the code-tracking (DLL) range error standard deviation for a
  /// given carrier-to-noise ratio.
  ///
  /// Used by GNSSMeasurements::ComputeSigmaRange to derive the pseudorange
  /// noise (`sigma_pseudorange_m`) stored on a GnssChannel from its
  /// CN0, by scaling the chip-fraction error returned here by `C * Tc`.
  ///
  /// @param params DLL tracking-loop parameters (bandwidth, integration time,
  ///                chip duration, spreading-code period, front-end bandwidth)
  /// @param CN0_w  Carrier-to-noise ratio [W/W, linear, not dB-Hz]
  /// @return       DLL tracking error standard deviation [chips]
  Real SigmaDll(const DllParams& params, Real CN0_w);

  /// @brief Compute the carrier-tracking (PLL) phase error standard deviation
  /// for a given carrier-to-noise ratio.
  ///
  /// Used by GNSSMeasurements::ComputeSigmaCarrierPhase to derive the
  /// carrier-phase noise (`sigma_carrier_phase_cycles`) stored on a
  /// GnssChannel, by scaling the result by `lambda / (2*pi)`.
  ///
  /// @param params PLL tracking-loop parameters (bandwidth, integration time)
  /// @param CN0_w  Carrier-to-noise ratio [W/W, linear, not dB-Hz]
  /// @return       PLL tracking error standard deviation [rad]
  Real SigmaPll(const PllParams& params, Real CN0_w);

  /// @brief Compute the frequency-tracking (FLL) error standard deviation for
  /// a given carrier-to-noise ratio.
  ///
  /// Used together with the carrier wavelength by
  /// GNSSMeasurements::ComputeSigmaRangeRate to derive the Doppler/range-rate
  /// noise (`sigma_doppler_hz`) stored on a GnssChannel.
  ///
  /// @param params FLL tracking-loop parameters (bandwidth, integration time,
  ///                CN0 threshold below which the squaring-loss factor doubles)
  /// @param CN0_w  Carrier-to-noise ratio [W/W, linear, not dB-Hz]
  /// @return       FLL tracking error standard deviation [Hz]
  Real SigmaFll(const FllParams& params, Real CN0_w);

  /// @brief Array (element-wise) overload of SigmaDll() over a vector of CN0 values.
  ArrX SigmaDll(const DllParams& params, const ArrX& CN0_w);
  /// @brief Array (element-wise) overload of SigmaPll() over a vector of CN0 values.
  ArrX SigmaPll(const PllParams& params, const ArrX& CN0_w);
  /// @brief Array (element-wise) overload of SigmaFll() over a vector of CN0 values.
  ArrX SigmaFll(const FllParams& params, const ArrX& CN0_w);

  /// @brief Compute the free-space path loss for a link of a given range and
  /// signal frequency.
  ///
  /// Used by the internal `LinkBudget` helper in gnss_measurement.cc (called
  /// from GNSSMeasurements::ComputeCN0) to convert the transmitter-to-receiver
  /// range into a loss term subtracted from the link budget when estimating CN0.
  ///
  /// @param dist Slant range between transmitter and receiver [m]
  /// @param freq Signal frequency [Hz]
  /// @return     Free-space path loss [dB]
  Real FreeSpacePathLoss(Real dist, Real freq);

  /// @brief Array (element-wise) overload of FreeSpacePathLoss() over a vector of ranges.
  ArrX FreeSpacePathLoss(const ArrX& dist, Real freq);

  /// @brief Evaluate a parabolic-dish antenna gain pattern at an off-boresight
  /// angle, using the standard parabolic approximation `G_max - 12*(phi/hpbw)^2`.
  ///
  /// Provides a simple analytic antenna-gain model (alternative to the
  /// measured patterns loaded by Antenna::ComputeGain) for link-budget /
  /// CN0-style calculations involving a parabolic-dish antenna.
  ///
  /// @param phi  Off-boresight (pointing) angle [same angular unit as `hpbw`,
  ///              e.g. deg or rad -- the two must be consistent]
  /// @param hpbw Half-power beamwidth [same angular unit as `phi`]
  /// @param G_max Peak (boresight) antenna gain [dB]
  /// @return     Antenna gain at `phi` [dB]
  Real ParabolicAntennaGain(Real phi, Real hpbw, Real G_max);

  /// @brief Array (element-wise) overload of ParabolicAntennaGain() over a vector of angles.
  ArrX ParabolicAntennaGain(const ArrX& phi, Real hpbw, Real G_max);

  /// @brief Line-of-sight visibility test between two points in a common
  /// frame, accounting for a spherical occluding body.
  ///
  /// Used by BuildChannels() (when `options.apply_visibility` is true) to
  /// drop transmitter channels whose line of sight to the receiver is
  /// blocked by an occluding body (e.g. the Moon or Earth), or whose
  /// elevation from a surface point falls below a minimum elevation mask.
  /// Handles three cases: either endpoint near the occluding body's surface
  /// (elevation-mask test) or both endpoints elevated (geometric horizon
  /// occlusion test).
  ///
  /// @param r1     First point [m], any common frame
  /// @param r2     Second point [m], same frame as `r1`
  /// @param R_body Radius of the occluding body [m]
  /// @param r_body Center of the occluding body [m], same frame as `r1`
  /// @param min_alt Minimum altitude for visibility test [m]
  /// @param min_elev_rad Minimum elevation angle for visibility test [rad]
  /// @return       True if `r1` and `r2` have an unobstructed line of sight
  ///                (subject to the elevation mask)
  bool ComputeVisibility(const Vec3& r1, const Vec3& r2, Real R_body,
                         const Vec3& r_body = Vec3::Zero(), const Real min_alt = 10e3,
                         const Real min_elev_deg = 5.0);

}  // namespace lupnt
