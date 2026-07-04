#include "lupnt/simulations/Ephemeris/ephemeris_simulation.h"

#include <algorithm>
#include <cmath>

#include "lupnt/applications/almanac.h"
#include "lupnt/applications/ephemeris.h"
#include "lupnt/lupnt.h"

namespace lupnt {

  namespace {

    Ptr<NBodyDynamics> MakeTruthDynamics(const EphemerisSimulationConfig& cfg) {
      auto dynamics = MakePtr<NBodyDynamics>();
      dynamics->SetIntegrator(IntegratorType::RKF45);
      // NOTE: tolerances tighter than ~1e-11 have been observed to make RKF45's
      // adaptive step-size control silently under-converge for highly eccentric
      // orbits (e.g. e=0.6 near periapsis), producing a plausible-looking but
      // wrong trajectory without throwing. These match the tolerances used
      // successfully in python/examples/ex1_propagate_orbit.ipynb.
      dynamics->SetIntegratorParams(IntegratorParams(20, 1.0e-6, 1.0e-10));
      dynamics->AddBody(Body::Moon(cfg.moon_gravity_degree, cfg.moon_gravity_order));
      if (cfg.include_earth) dynamics->AddBody(Body::Earth());
      if (cfg.include_sun) dynamics->AddBody(Body::Sun());
      dynamics->SetFrame(cfg.propagate_frame);
      dynamics->SetTimeStep(cfg.integration_step_s);
      dynamics->SetUseRelativity(cfg.use_relativity);
      return dynamics;
    }

    // Number of bits needed to represent a parameter with the given (non-negative)
    // `range` at a resolution of 2^-bits_frac, plus 1 sign bit if `signed_range`
    // and 1 margin bit if `margin_bit` (mirrors common GNSS broadcast-message
    // bit-allocation conventions, e.g. GPS LNAV/CNAV scale-factor tables).
    int RangeBits(double range, int bits_frac, bool signed_range, bool margin_bit) {
      const int bits_range
          = range > 0.0 ? std::max(1, static_cast<int>(std::ceil(std::log2(range)))) : 1;
      return bits_range + std::max(bits_frac, 0) + (signed_range ? 1 : 0) + (margin_bit ? 1 : 0);
    }

    // Binary search for the minimum number of fractional bits `k` such that
    // quantizing parameter `idx` of `params` to a resolution of 2^-k changes the
    // fitted trajectory (evaluated at `t_s`, mean position offset) by less than
    // `precision_m`, holding every other parameter exact. Mirrors the per-
    // parameter resolution search used to size GNSS broadcast ephemeris fields.
    template <typename EphemT>
    int SearchParamBits(const EphemT& eph, const VecXd& t_s, const VecXd& params, int idx,
                        double precision_m, int max_bits = 50) {
      const MatXd baseline = eph.Eval(t_s, params);
      auto perturbed_error = [&](int k) {
        const double eps = std::pow(2.0, -k);
        VecXd p = params;
        p(idx) += eps;
        const MatXd plus = eph.Eval(t_s, p);
        p(idx) -= 2.0 * eps;
        const MatXd minus = eph.Eval(t_s, p);
        const double err_plus = (plus.leftCols(3) - baseline.leftCols(3)).rowwise().norm().mean();
        const double err_minus = (minus.leftCols(3) - baseline.leftCols(3)).rowwise().norm().mean();
        return std::max(err_plus, err_minus);
      };
      int lo = 0, hi = max_bits;
      if (perturbed_error(lo) <= precision_m) return lo;
      if (perturbed_error(hi) > precision_m) return hi;
      while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (perturbed_error(mid) <= precision_m) {
          hi = mid;
        } else {
          lo = mid;
        }
      }
      return hi;
    }

    // True if parameter `name` should be treated as a 2*pi-ranged, signed angle
    // for bit-allocation purposes (both the raw Keplerian angles used by
    // CartesianEphemeris, and the Chebyshev coefficients of the corresponding
    // time-varying elements used by Almanac -- a conservative approximation for
    // the higher-order coefficients, which are typically small in magnitude).
    bool IsAngleParam(const std::string& name) {
      return name.rfind("i", 0) == 0 || name.rfind("raan", 0) == 0 || name.rfind("argp", 0) == 0
             || name.rfind("M", 0) == 0;
    }
    bool IsEccentricityParam(const std::string& name) { return name.rfind("e", 0) == 0; }
    bool IsSemiMajorAxisParam(const std::string& name) { return name.rfind("a", 0) == 0; }

    template <typename EphemT>
    EphemerisWindowResult SummarizeWindows(const EphemT& eph, const std::vector<VecXd>& t_windows,
                                            const std::vector<MatXd>& rv_windows,
                                            double fit_window_min, double precision_m) {
      const int num_params = eph.NumParams();
      const std::vector<std::string> names = eph.ParamNames();
      VecXd max_k = VecXd::Zero(num_params);
      VecXd max_abs = VecXd::Zero(num_params);

      double pos_rms_sum = 0.0, vel_rms_sum = 0.0, pos_p95_sum = 0.0, vel_p95_sum = 0.0;

      for (size_t w = 0; w < t_windows.size(); ++w) {
        const VecXd params = eph.Fit(t_windows[w], rv_windows[w]);
        const EphemerisFitErrorStats stats = eph.EvalError(t_windows[w], rv_windows[w], params);
        pos_rms_sum += stats.rms_pos_m(3);
        vel_rms_sum += stats.rms_vel_mps(3);
        pos_p95_sum += stats.p95_pos_m(3);
        vel_p95_sum += stats.p95_vel_mps(3);

        for (int j = 2; j < num_params; ++j) {
          max_abs(j) = std::max(max_abs(j), std::abs(params(j)));
          const int k = SearchParamBits(eph, t_windows[w], params, j, precision_m);
          max_k(j) = std::max(max_k(j), static_cast<double>(k));
        }
      }

      const int nw = static_cast<int>(t_windows.size());
      EphemerisWindowResult result;
      result.fit_window_min = fit_window_min;
      result.num_params = num_params;
      result.pos_rms_m = pos_rms_sum / nw;
      result.vel_rms_mps = vel_rms_sum / nw;
      result.pos_p95_m = pos_p95_sum / nw;
      result.vel_p95_mps = vel_p95_sum / nw;

      int total_bits = 0;
      for (int j = 2; j < num_params; ++j) {
        double range;
        bool signed_range, margin_bit;
        if (IsAngleParam(names[j])) {
          range = 2.0 * PI;
          signed_range = true;
          margin_bit = false;
        } else if (IsEccentricityParam(names[j])) {
          range = 1.0;
          signed_range = false;
          margin_bit = false;
        } else if (IsSemiMajorAxisParam(names[j])) {
          range = max_abs(j);
          signed_range = false;
          margin_bit = true;
        } else {
          range = max_abs(j);
          signed_range = true;
          margin_bit = true;
        }
        total_bits += RangeBits(range, static_cast<int>(max_k(j)), signed_range, margin_bit);
      }
      result.total_bits = total_bits;
      return result;
    }

  }  // namespace

  EphemerisSimulation::EphemerisSimulation(EphemerisSimulationConfig config)
      : config_(std::move(config)) {}

  void EphemerisSimulation::Setup() {
    LUPNT_CHECK(config_.duration_days > 0.0, "duration_days must be positive", "EphemerisSimulation");
    LUPNT_CHECK(config_.sample_dt_s > 0.0, "sample_dt_s must be positive", "EphemerisSimulation");
    LUPNT_CHECK(!config_.fit_window_minutes.empty(), "fit_window_minutes must not be empty",
                "EphemerisSimulation");
    LUPNT_CHECK(config_.num_windows >= 2, "num_windows must be at least 2", "EphemerisSimulation");
    LUPNT_CHECK(config_.datasize_precision_m > 0.0, "datasize_precision_m must be positive",
                "EphemerisSimulation");
    setup_complete_ = true;
  }

  void EphemerisSimulation::Run() {
    LUPNT_CHECK(setup_complete_, "Call Setup() before Run()", "EphemerisSimulation");
    const EphemerisSimulationConfig& cfg = config_;

    const Real t0_tdb = ConvertTime(GregorianToTime(cfg.start_epoch_utc), Time::UTC, Time::TDB);

    Vec6 coe0;
    coe0 << cfg.orbit.a_m, cfg.orbit.ecc, cfg.orbit.inc_rad, cfg.orbit.raan_rad, cfg.orbit.argp_rad,
        cfg.orbit.m0_rad;
    const Vec6 rv0_coe_frame = ClassicalToCart(coe0, Real(GM_MOON));
    const Vec6 rv0 = cfg.orbit.coe_frame == cfg.propagate_frame
                         ? rv0_coe_frame
                         : ConvertFrame(t0_tdb, rv0_coe_frame, cfg.orbit.coe_frame,
                                        cfg.propagate_frame);

    Ptr<NBodyDynamics> dynamics = MakeTruthDynamics(cfg);

    const int n
        = static_cast<int>(std::floor(cfg.duration_days * SECS_DAY / cfg.sample_dt_s + 1.0e-9)) + 1;
    VecX tspan(n);
    for (int i = 0; i < n; ++i) tspan(i) = i * cfg.sample_dt_s;
    VecX tfs(n);
    for (int i = 0; i < n; ++i) tfs(i) = t0_tdb + tspan(i);

    Cart6 x0(rv0, cfg.propagate_frame);
    const MatX rv_prop = dynamics->Propagate(x0, tfs);

    t_truth_s_ = tspan.cast<double>();
    rv_truth_ = rv_prop.cast<double>();

    cartesian_results_.clear();
    almanac_results_.clear();

    const CartesianEphemeris cart_eph(EphemerisFitOptions{
        cfg.cartesian_poly_order, cfg.cartesian_use_keplerian_baseline, GM_MOON});
    const Almanac almanac(AlmanacFitOptions{cfg.almanac_poly_order, GM_MOON});

    auto pbar = Logger::GetProgressBar(static_cast<int>(cfg.fit_window_minutes.size()),
                                       "Ephemeris datasize/accuracy sweep");
    int prog = 0;
    for (double fit_min : cfg.fit_window_minutes) {
      const double window_s = fit_min * SECS_MINUTE;
      const int window_points = static_cast<int>(std::round(window_s / cfg.sample_dt_s)) + 1;
      LUPNT_CHECK(window_points >= 2 && window_points <= n,
                  "fit_window_minutes must fit within the sampled truth trajectory",
                  "EphemerisSimulation");

      const int margin = window_points / 2;
      const int span = std::max(1, n - 1 - 2 * margin);
      std::vector<VecXd> t_windows;
      std::vector<MatXd> rv_windows;
      for (int w = 0; w < cfg.num_windows; ++w) {
        const double frac
            = cfg.num_windows > 1 ? static_cast<double>(w) / (cfg.num_windows - 1) : 0.5;
        int start = margin + static_cast<int>(std::round(frac * span)) - window_points / 2;
        start = std::max(0, std::min(start, n - window_points));
        t_windows.push_back(t_truth_s_.segment(start, window_points));
        rv_windows.push_back(rv_truth_.middleRows(start, window_points));
      }

      cartesian_results_.push_back(
          SummarizeWindows(cart_eph, t_windows, rv_windows, fit_min, cfg.datasize_precision_m));
      almanac_results_.push_back(
          SummarizeWindows(almanac, t_windows, rv_windows, fit_min, cfg.datasize_precision_m));

      pbar->Update(++prog);
    }
    pbar->Finish();
  }

}  // namespace lupnt
