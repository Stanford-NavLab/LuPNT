// Example 2: Relativistic Time Conversions -- TT -> TCL -> LT
// ------------------------------------------------------------------------------
// C++ counterpart of python/examples/ex2_time_conversions.ipynb (Section 4.1).
//
// General relativity gives every gravitational environment its own proper time.
// The IAU time hierarchy relevant here is:
//
//   TT   Terrestrial Time      (geocentric; TAI + 32.184 s realises it)
//   TCG  Geocentric Coord Time  (GCRS coordinate time)
//   TCB  Barycentric Coord Time (BCRS coordinate time)
//   TDB  Barycentric Dynamical Time
//   TCL  Selenocentric Coord Time
//   LT   Lunar Time            (proper time on the lunar geoid; a.k.a. TL)
//
// The defining IAU scale factors (all in lupnt/core/constants.h) are:
//   L_G  TCG->TT rate, L_B TCB->TDB rate, L_L TCL->TL rate,
//   L_EM mean LT-TT rate (Fienga et al. 2023).
//
// This example reproduces the analytical *secular drift rates* between the
// scales; the TCG-TT and TL-TCL relations are exact-linear (no periodic terms).

#include <iomanip>
#include <iostream>

#include "lupnt/lupnt.h"

using namespace lupnt;

int main() {
  std::cout << std::scientific << std::setprecision(6);
  std::cout << "IAU / lunar scale factors:\n";
  std::cout << "  L_G  = " << L_G << "  (TCG->TT)\n";
  std::cout << "  L_B  = " << L_B << "  (TCB->TDB)\n";
  std::cout << "  L_L  = " << L_L << "  (TCL->TL)\n";
  std::cout << "  L_EM = " << L_EM << "  (mean LT-TT)\n";
  std::cout << "  c    = " << C << " m/s\n\n";

  // All scales coincide at the reference epoch T_0 = 1977-01-01 00:00:32.184 TAI.
  std::cout << "Reference epoch T_0: MJD " << MJD_COORDINATE_TT_TCG_TCB << " (TAI)\n\n";

  // --- Section 4.1: analytical secular drift rates [microseconds / day] ------
  const double us_per_day = SECS_DAY * 1e6;

  // TCG - TT: exact linear (IAU Resolution B1.5), d(TCG-TT)/dTT = L_G/(1-L_G).
  const double rate_tcg_tt = L_G / (1.0 - L_G) * us_per_day;

  // TL - TCL: exact linear on the lunar geoid, d(TL-TCL)/dTCL = -L_L.
  const double rate_lt_tcl = -L_L * us_per_day;

  // LT - TT secular rate: combines the TCL-TT and LT-TCL rates with the mean
  // tidal correction L_EM (Fienga et al. 2023).
  const double rate_lt_tt = ((L_G - L_L) / (1.0 - L_B) - L_EM) * us_per_day;

  std::cout << std::fixed << std::setprecision(4);
  std::cout << "Secular drift rates:\n";
  std::cout << "  d(TCG - TT)/dTT   = " << rate_tcg_tt << " us/day\n";
  std::cout << "  d(TL  - TCL)/dTCL = " << rate_lt_tcl << " us/day\n";
  std::cout << "  d(LT  - TT)/dTT   = " << rate_lt_tt << " us/day\n";

  // Accumulated LT-TT offset over a one-year mission arc.
  std::cout << "\nAccumulated LT-TT drift over 365.25 days: " << rate_lt_tt * 365.25 << " us\n";

  return 0;
}
