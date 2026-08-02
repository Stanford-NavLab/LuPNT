// Example 2: Relativistic Time Conversions in LuPNT
// -----------------------------------------------------------------------------
// C++ counterpart of python/examples/ex2_time_conversions.ipynb.
//
// HOW TIME CONVERSION WORKS IN LUPNT
//
// A time scale is a coordinate on a particular 4-D reference system, so
// converting between two of them means evaluating a relativistic relationship,
// not adding a constant. LuPNT organises this as a GRAPH: each scale is a node,
// and each registered edge returns the OFFSET between two scales,
//
//     delta_{A->B}(t) = t_B - t_A
//
// given the epoch's reading in scale A. A conversion walks the shortest
// registered route and sums the offsets along it.
//
//     constant   TAI<->TT (32.184 s), TAI<->GPS (19 s)
//     table      TAI<->UTC (leap seconds), UTC<->UT1 (Earth orientation)
//     linear     TT<->TCG, TDB<->TCB, TCL<->LT   (rescale by L ~ 1e-8)
//     model      TT<->TDB                        (Chebyshev fit / series)
//     integral   TDB<->TCL                       (sweep from T_0 = 1977)
//
// Two consequences matter in practice.
//
// OFFSETS, NOT ABSOLUTE EPOCHS. An absolute epoch is a double counting seconds
// from J2000, so near 2030 one unit in the last place is
//
//     |t| * 2^-52 ~ 2.45e-7 s ~ 245 ns ~ 73 m * c.
//
// Summing small offsets never differences two large numbers, so Epoch keeps full
// double precision. ConvertTime has to RETURN an absolute epoch and is therefore
// capped near that floor however good the model is. Use Epoch, or the offset
// accessors (TtMinusTdb / TdbMinusTcl / TdbMinusLt), whenever the difference
// itself is the quantity of interest.
//
// ROUTING MATTERS. TCL<->LT is a pure linear rescaling; reaching it via TAI
// would drag in the TDB<->TCL integral, which sweeps from 1977 on every call.
// The graph uses the direct edge instead.
//
// This example compares the three TT-TDB routes LuPNT offers, then uses the
// integral-based TDB-LT relation to reach TT, and finally decomposes the
// periodic structure.

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

#include "lupnt/lupnt.h"

using namespace lupnt;

namespace {

  constexpr double kSecsYear = 365.25 * 86400.0;

  /// Constant + linear drift + what is left.
  ///
  /// Every residual here turns out to be almost purely constant+linear, so the
  /// drift and the DETRENDED rms are the informative numbers: a raw rms mostly
  /// reports the constant offset and hides whether the physics is right.
  struct Decomposition {
    double rms, constant, drift, detrended_rms;
  };

  Decomposition Decompose(const std::vector<double>& d, const std::vector<double>& t) {
    const int n = static_cast<int>(d.size());
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (int i = 0; i < n; i++) {
      const double x = (t[i] - t[0]) / kSecsYear;
      sx += x;
      sy += d[i];
      sxx += x * x;
      sxy += x * d[i];
    }
    const double den = n * sxx - sx * sx;
    const double slope = (n * sxy - sx * sy) / den;
    const double icept = (sy * sxx - sx * sxy) / den;
    double s2 = 0, r2 = 0;
    for (int i = 0; i < n; i++) {
      const double x = (t[i] - t[0]) / kSecsYear;
      s2 += d[i] * d[i];
      const double r = d[i] - (icept + slope * x);
      r2 += r * r;
    }
    return {std::sqrt(s2 / n), icept, slope, std::sqrt(r2 / n)};
  }

  /// Amplitude spectrum of a detrended, uniformly-sampled signal.
  /// Naive DFT -- the sample counts here do not justify pulling in an FFT.
  void Spectrum(const std::vector<double>& sig, double dt_s, std::vector<double>* period_d,
                std::vector<double>* amp, std::vector<double>* freq) {
    const int n = static_cast<int>(sig.size());
    double mean = 0.0;
    for (double v : sig) mean += v;
    mean /= n;

    std::vector<double> x(n);
    double wsum = 0.0;
    for (int i = 0; i < n; i++) {
      const double w = 0.5 * (1.0 - std::cos(2.0 * PI * i / (n - 1)));  // Hanning
      x[i] = (sig[i] - mean) * w;
      wsum += w;
    }
    period_d->clear();
    amp->clear();
    freq->clear();
    for (int k = 1; k < n / 2; k++) {
      std::complex<double> acc(0.0, 0.0);
      for (int i = 0; i < n; i++) {
        const double ph = -2.0 * PI * k * i / n;
        acc += x[i] * std::complex<double>(std::cos(ph), std::sin(ph));
      }
      const double f = static_cast<double>(k) / (n * dt_s);  // [1/s]
      freq->push_back(f);
      period_d->push_back(1.0 / (f * SECS_DAY));
      amp->push_back(std::abs(acc) / (wsum / 2.0));
    }
  }

  /// Strongest lines in a period band, one per resolution element.
  ///
  /// Spectral leakage puts sidelobes about ONE frequency bin either side of a
  /// real line, so the separation floor has to be expressed in BINS -- a fixed
  /// log-period spacing is too wide at short periods and too narrow at long
  /// ones. A Hanning mainlobe spans ~4 bins; 3 rejects its skirts while keeping
  /// genuinely distinct lines.
  std::vector<std::pair<double, double>> Peaks(const std::vector<double>& period_d,
                                               const std::vector<double>& amp,
                                               const std::vector<double>& freq, double pmin,
                                               double pmax, int want, double sep_bins = 3.0) {
    const double df = freq[1] - freq[0];
    std::vector<int> idx;
    for (size_t i = 0; i < period_d.size(); i++)
      if (period_d[i] >= pmin && period_d[i] <= pmax) idx.push_back(static_cast<int>(i));
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return amp[a] > amp[b]; });

    std::vector<std::pair<double, double>> out;
    std::vector<double> taken;
    for (int i : idx) {
      bool clear = true;
      for (double f : taken)
        if (std::abs(freq[i] - f) <= sep_bins * df) clear = false;
      if (!clear) continue;
      out.emplace_back(period_d[i], amp[i]);
      taken.push_back(freq[i]);
      if (static_cast<int>(out.size()) == want) break;
    }
    return out;
  }

}  // namespace

int main() {
  // Build the grid EXACTLY: Epoch::FromGregorian uses integer arithmetic, so the
  // instant carries no input rounding. GregorianToTime would route via a
  // Modified Julian Date, whose 1858 origin makes its ULP ~1.1 us at present-day
  // epochs -- 8x coarser than seconds-from-J2000.
  const Epoch e_start = Epoch::FromGregorian(2020, 1, 1, 0, 0, Real(0.0), Time::TDB);
  const Epoch e_end = Epoch::FromGregorian(2035, 1, 1, 0, 0, Real(0.0), Time::TDB);
  const int n_pts = 721;  // ~7.6 day sampling over 15 years

  const double t0 = e_start.ToSeconds().val();
  const double t1 = e_end.ToSeconds().val();
  std::vector<double> t(n_pts);
  for (int i = 0; i < n_pts; i++) t[i] = t0 + (t1 - t0) * i / (n_pts - 1);

  std::cout << std::fixed << std::setprecision(1);
  std::cout << "grid: " << n_pts << " points, 2020 - 2035\n"
            << "absolute-epoch ULP here: " << std::abs(t1) * std::pow(2.0, -52) * 1e9
            << " ns  (every route below returns OFFSETS)\n\n";

  // ---------------------------------------------------------------------------
  // 1. TT - TDB: three routes
  // ---------------------------------------------------------------------------
  // TDB is the barycentric coordinate time JPL publishes its ephemerides in; TT
  // is the terrestrial scale clocks on the geoid realise. Their difference is
  // dominated by a ~1.66 ms annual term driven by the Earth's eccentric orbit:
  // the Earth's speed and its depth in the Sun's potential both vary over the
  // year, so a terrestrial clock gains and loses against a barycentric one.
  std::vector<double> tt_tdb_kernel(n_pts), tt_tdb_integral(n_pts), tt_tdb_fit(n_pts);

  // Route 2 (the reference): JPL's own integrated TT-TDB, read from the de440t
  // time-ephemeris segment as a clean offset.
  for (int i = 0; i < n_pts; i++)
    tt_tdb_kernel[i] = spice::GetTimeEphemerisOffset(Real(t[i]), spice::kNaifTtMinusTdb).val();

  // Route 1: the DE440 Eq. (3) relativistic integral -- the definition,
  // integrated from LuPNT's own ephemeris rather than reading JPL's answer.
  SetDe440TtTdbStep(0.01 * SECS_DAY);  // 864 s trapezoid
  {
    VecX tv(n_pts);
    for (int i = 0; i < n_pts; i++) tv(i) = t[i];
    const VecX out = TdbMinusTtDe440(tv);  // single sorted sweep, not N integrations
    for (int i = 0; i < n_pts; i++) tt_tdb_integral[i] = -out(i).val();
  }

  // Route 3: the Chebyshev fit of that kernel -- what LuPNT uses by default.
  InitTtMinusTdbFit(Real(t0 - 30 * SECS_DAY), Real(t1 + 30 * SECS_DAY));
  for (int i = 0; i < n_pts; i++) tt_tdb_fit[i] = TtMinusTdb(Real(t[i])).val();

  std::cout << std::left << std::setw(34) << "route" << std::right << std::setw(12) << "rms"
            << std::setw(15) << "drift" << std::setw(16) << "detrended rms" << "\n"
            << std::string(77, '-') << "\n";
  const std::pair<const char*, const std::vector<double>*> routes[] = {
      {"DE440 Eq.(3) integral", &tt_tdb_integral},
      {"Chebyshev fit of the kernel", &tt_tdb_fit},
  };
  for (const auto& [name, v] : routes) {
    std::vector<double> r(n_pts);
    for (int i = 0; i < n_pts; i++) r[i] = ((*v)[i] - tt_tdb_kernel[i]) * 1e9;  // [ns]
    const Decomposition d = Decompose(r, t);
    std::cout << std::left << std::setw(34) << name << std::right << std::fixed
              << std::setprecision(4) << std::setw(9) << d.rms << " ns" << std::showpos
              << std::setw(11) << d.drift << std::noshowpos << " ns/yr" << std::setw(12)
              << d.detrended_rms << " ns\n";
  }
  std::cout << "\n  The fit reproduces the kernel to sub-picosecond -- it interpolates the same\n"
               "  data. Eq.(3) is an independent computation, and its difference is almost\n"
               "  purely SECULAR: the periodic physics matches, only the rate differs. That\n"
               "  residual rate comes from the small bodies DE440 integrates (343 asteroids,\n"
               "  30 KBOs and a 36-point Kuiper ring at 44 au), modelled here as rings using\n"
               "  DE440's own header constants.\n\n";

  // ---------------------------------------------------------------------------
  // 2. TDB - LT by integration, then on to TT
  // ---------------------------------------------------------------------------
  // LT is the analogue of TT for a clock on the lunar selenoid. Reaching it from
  // TDB uses the one genuinely expensive edge in the graph -- the TDB<->TCL
  // integral -- followed by the cheap linear TCL<->LT rescaling by L_L:
  //
  //     LT - TT = (LT - TDB) + (TDB - TT)
  //
  // which is what Epoch composes when asked for LT -> TT. Every term is a small
  // offset, so the ~245 ns epoch floor never enters.
  std::vector<double> lt_tt(n_pts);
  for (int i = 0; i < n_pts; i++) lt_tt[i] = -TdbMinusLt(Real(t[i])).val() - tt_tdb_fit[i];

  // Cross-check against Epoch walking the graph itself. Mind the sign:
  // TimeScaleOffset(Epoch(LT), TT) returns TT - LT.
  double max_diff = 0.0;
  for (int i = 0; i < n_pts; i += std::max(1, n_pts / 40)) {
    const Epoch e_lt = Epoch::FromSeconds(Real(t[i] - TdbMinusLt(Real(t[i])).val()), Time::LT);
    const double via_epoch = -TimeScaleOffset(e_lt, Time::TT).val();
    max_diff = std::max(max_diff, std::abs(lt_tt[i] - via_epoch));
  }
  std::cout << "chain vs Epoch(LT->TT): max |diff| = " << std::setprecision(4) << max_diff * 1e9
            << " ns\n";
  {
    const Decomposition d = Decompose(lt_tt, t);
    std::cout << "secular drift of LT - TT: " << std::setprecision(4) << d.drift * 1e6 / 365.25
              << " us/day\n  (Turyshev et al. 2025 give 56.0256 us/day)\n\n";
  }

  // ---------------------------------------------------------------------------
  // 3. Fourier decomposition
  // ---------------------------------------------------------------------------
  // With the secular trend removed, what remains is a sum of periodic terms
  // whose frequencies are the orbital periods driving each effect.
  //
  // SAMPLING GOVERNS WHAT CAN BE SEEN. The 15-year grid has ~7.6 day spacing, so
  // its Nyquist period is ~15 days -- ample for the Earth's annual term (~48
  // samples per period), but the lunar terms at 27.55 d (anomalistic month) and
  // 29.53 d (synodic month) get under 4 samples each and alias badly.
  //
  // So the two signals need different grids. Resolution is set by record LENGTH,
  // aliasing by sample SPACING: the two lunar months sit 2.4e-3 /d apart in
  // frequency, so a record of T days separates them by 2.4e-3*T bins. A Hanning
  // mainlobe spans ~4 bins, so T must exceed ~4.5 years -- 3 years is not enough
  // and merges them however the peaks are picked.
  {
    const Decomposition d = Decompose(tt_tdb_kernel, t);
    std::vector<double> periodic(n_pts);
    for (int i = 0; i < n_pts; i++) {
      const double x = (t[i] - t[0]) / kSecsYear;
      periodic[i] = tt_tdb_kernel[i] - (d.constant + d.drift * x);
    }
    std::vector<double> p, a, f;
    Spectrum(periodic, t[1] - t[0], &p, &a, &f);
    std::cout << "TT - TDB   (" << n_pts << " pts over 15 yr, " << std::setprecision(1)
              << (t[1] - t[0]) / SECS_DAY << " d spacing)\n";
    // Cap the band at ~1/3 of the record: longer "periods" are the residual
    // trend leaking, not planetary terms.
    for (const auto& [per, amp] : Peaks(p, a, f, 100.0, 5.0 * 365.25, 4))
      std::cout << "   " << std::setw(8) << std::setprecision(3) << per / 365.25
                << " yr   amplitude " << std::setw(9) << std::setprecision(2) << amp * 1e6
                << " us\n";
  }

  {
    const double years = 6.0;
    const int n_lun = static_cast<int>(years * 365.25 * 4) + 1;  // 0.25 day spacing
    std::vector<double> t_lun(n_lun), lt(n_lun);
    for (int i = 0; i < n_lun; i++) t_lun[i] = t0 + years * 365.25 * SECS_DAY * i / (n_lun - 1);
    InitLtMinusTtFit(Real(t_lun[0] - 30 * SECS_DAY), Real(t_lun[n_lun - 1] + 30 * SECS_DAY));
    for (int i = 0; i < n_lun; i++)
      lt[i] = -TdbMinusLt(Real(t_lun[i])).val() - TtMinusTdb(Real(t_lun[i])).val();

    const Decomposition d = Decompose(lt, t_lun);
    std::vector<double> periodic(n_lun);
    for (int i = 0; i < n_lun; i++) {
      const double x = (t_lun[i] - t_lun[0]) / kSecsYear;
      periodic[i] = lt[i] - (d.constant + d.drift * x);
    }
    std::vector<double> p, a, f;
    Spectrum(periodic, t_lun[1] - t_lun[0], &p, &a, &f);
    const double sep = std::abs(1.0 / 29.5306 - 1.0 / 27.5545) * years * 365.25;
    std::cout << "\nLT - TT    (" << n_lun << " pts over " << std::setprecision(0) << years
              << " yr, 0.25 d spacing; the two lunar months are " << std::setprecision(1) << sep
              << " bins apart)\n";
    for (const auto& [per, amp] : Peaks(p, a, f, 5.0, 200.0, 5))
      std::cout << "   " << std::setw(8) << std::setprecision(3) << per << " d    amplitude "
                << std::setw(9) << std::setprecision(4) << amp * 1e6 << " us\n";
  }

  std::cout << "\nChoosing a route:\n"
               "  sub-ns, or a time multiplied by c    -> Epoch, or the offset accessors\n"
               "  agreement with JPL to picoseconds    -> the Chebyshev fit (the default)\n"
               "  the definition from first principles -> SetTtTdbModel(DE440_INTEGRAL)\n"
               "  microsecond tolerances, convenience  -> ConvertTime (245 ns epoch floor)\n";
  return 0;
}
