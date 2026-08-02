/**
 * @file gcpm/gcpm_surrogate.h
 * @brief Fast interpolation surrogate for GCPM electron density.
 *
 * GCPM costs ~10 IRI evaluations per call (~0.5 ms). This surrogate replaces the
 * whole model with a 7-D interpolation of pre-sampled GCPM output in its natural
 * solar-magnetic frame plus UT:
 *
 *     Ne = f(r, lam_m, MLT, UT, Kp, R12, doy)
 *
 * The table (log Ne, float32) is generated offline by
 * scripts/gcpm_surrogate_sample.py + gcpm_surrogate_build.py and loaded at
 * runtime from a binary file. Interpolation is multilinear in log-density with
 * MLT/UT/doy treated as periodic. Reproduces GCPM to <1% in the smooth
 * plasmasphere/topside and a few % near the sharp F2 peak, at >10x the speed.
 *
 * @copyright Copyright (c) 2026
 */
#pragma once

#include <string>
#include <vector>

namespace pecsim {

  class GcpmSurrogate {
  public:
    /** Load the surrogate table from a binary file. Returns false if the file is
     *  absent/unreadable (caller should then fall back to full GCPM). */
    bool load(const std::string& path);

    bool loaded() const { return loaded_; }

    /** Interpolated GCPM electron density [cm^-3].
     *  @param r     geocentric radius [RE]
     *  @param lam_m SM magnetic latitude [deg]
     *  @param mlt   magnetic local time [h]
     *  @param ut    universal time [h]
     *  @param kp    geomagnetic index
     *  @param r12   sunspot number
     *  @param doy   day of year */
    double eval(double r, double lam_m, double mlt, double ut, double kp, double r12,
                double doy) const;

    /** Process-wide surrogate, lazily loaded from LUPNT_GCPM_SURROGATE or
     *  <plasma base>/gcpm_surrogate.bin on first use. */
    static const GcpmSurrogate& instance();

  private:
    bool loaded_ = false;
    std::vector<std::vector<double>> axes_;  // 7 axes: r, lamm, mlt, ut, kp, r12, doy
    std::vector<bool> periodic_;             // per-axis periodicity
    std::vector<double> period_;             // period (if periodic)
    std::vector<float> table_;               // flattened log(Ne), C-order over axes_
    std::vector<long> stride_;               // per-axis flat stride
  };

}  // namespace pecsim
