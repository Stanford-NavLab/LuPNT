/**
 * @file antenna.h
 * @author Stanford NAV LAB
 * @brief Antenna class
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <string>

#include "lupnt/core/definitions.h"
#include "lupnt/numerics/vector_macros.h"

namespace lupnt {
  class Antenna {
  public:
    Antenna() = default;
    Antenna(const std::string& name) : name_(name) { LoadAntennaPattern(); };

    /// @brief Load the antenna gain pattern from the named pattern file (e.g.
    /// "Block-IIR_ACE") into `gain_`/`phi_`/`theta_`, reformatting the raw
    /// 1D/2D table into LuPNT's [-180,180] (phi) x [0,360] (theta) convention.
    ///
    /// Called from the `Antenna(name)` constructor so that
    /// `GnssConstellation`/`GNSSMeasurements` (e.g.
    /// GNSSMeasurements::SetReceiverAntenna,
    /// GnssConstellation::GetTransmitterAntenna) get a ready-to-query gain
    /// pattern for CN0/link-budget calculations. If `name_` is empty, the
    /// antenna is treated as omni-directional (ComputeGain returns 0).
    void LoadAntennaPattern();

    /// @brief Look up the antenna gain at a given boresight-relative
    /// elevation/azimuth pair from the loaded pattern.
    ///
    /// Used by GNSSMeasurements::ComputeCN0 to evaluate the transmitter and
    /// receiver antenna gains (`G_tx`, `G_rx`) that feed the link-budget /
    /// CN0 estimate for a GNSS channel.
    ///
    /// @param theta Azimuth angle around boresight [rad], wrapped to [0, 2*pi)
    /// @param phi   Elevation/off-boresight angle [rad], wrapped to [-pi, pi]
    /// @return      Antenna gain at (`theta`, `phi`) [dB]; 0 for an
    ///              omni-directional antenna, NaN if `phi` exceeds the
    ///              pattern's coverage
    Real ComputeGain(Real theta, Real phi) const;
    VEC_DEF_REAL_REAL(ComputeGain)

    /// @brief Name of the loaded antenna pattern (e.g. "Block-IIR_ACE"); empty for
    /// omni-directional.
    std::string GetName() { return name_; }

    /// @brief Look up the antenna gain in the direction of a given unit
    /// vector, expressed in the antenna's own boresight-aligned frame.
    ///
    /// Convenience overload of ComputeGain(theta, phi) that derives the
    /// elevation/azimuth pair from a 3D `direction` vector (Z = boresight)
    /// before delegating to the 2-angle form.
    ///
    /// @param direction Direction to evaluate the gain in, in the antenna
    ///                   boresight frame (need not be normalized)
    /// @return Antenna gain in that direction [dB]
    Real ComputeGain(const Vec3& direction) const;

    /// @brief Raw 2D gain pattern table [dB], indexed as (phi, theta).
    MatXd GetGainMatrix() { return gain_; }

    /// @brief Phi (elevation/off-boresight) angle grid of the loaded pattern [rad].
    VecXd GetPhiVector() { return phi_; }

    /// @brief Theta (azimuth) angle grid of the loaded pattern [deg].
    VecXd GetThetaVector() { return theta_; }

  private:
    int n_dim_ = 0;          // Number of dimensions (0=omni, 1, or 2)
    std::string name_;       // Name (e.g., Block-IIR_ACE)
    double phi_max_ = 90.0;  // Maximum off-boresight angle [deg]
    MatXd gain_;             // Gain pattern [dB]
    VecXd phi_;              // Phi angles [deg]
    VecXd theta_;            // Theta angles [deg]

    /// @brief Normalize a raw 1D/2D antenna pattern table (from
    /// `LoadAntennaPattern`) into LuPNT's phi in [-180,180] deg,
    /// theta in [0,360] deg convention, padding/mirroring the gain table as
    /// needed so `ComputeGain` can interpolate over the full range.
    ///
    /// @param phi   Elevation/off-boresight angle samples [deg] (input grid;
    ///               rewritten in place to the normalized grid)
    /// @param theta Azimuth angle samples [deg] (input grid; rewritten in
    ///               place to the normalized [0,360] grid)
    /// @param gain  Gain table [dB], indexed as gain[phi][theta] (rewritten in
    ///               place to match the normalized `phi`/`theta` grids)
    void FormatAntennaPattern(std::vector<double>& phi, std::vector<double>& theta,
                              std::vector<std::vector<double>>& gain);
  };

}  // namespace lupnt
