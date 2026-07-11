#include "lupnt/applications/ephemeris/ephemeris_app.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "lupnt/agents/agent.h"
#include "lupnt/applications/ephemeris/lunanet_almanac.h"
#include "lupnt/applications/ephemeris/lunanet_ephemeris.h"
#include "lupnt/core/asset_factory.h"
#include "lupnt/lupnt.h"
#include "lupnt/simulations/simulation.h"

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
    // fitted trajectory (evaluated at `t_s`, *maximum* position offset over the
    // window) by less than `precision_m`, holding every other parameter exact.
    // Bounding the worst-case (rather than mean) offset keeps the realized
    // quantized error under control everywhere -- important for eccentric orbits,
    // where a coarse parameter step can spike the position near periapsis or
    // detune a period-dependent term. Mirrors the per-parameter resolution search
    // used to size GNSS broadcast ephemeris fields.
    template <typename EphemT> int SearchParamBits(const EphemT& eph, const VecXd& t_s,
                                                   const VecXd& params, int idx, double precision_m,
                                                   int max_bits = 50) {
      const MatXd baseline = eph.Eval(t_s, params);
      auto perturbed_error = [&](int k) {
        const double eps = std::pow(2.0, -k);
        VecXd p = params;
        p(idx) += eps;
        const MatXd plus = eph.Eval(t_s, p);
        p(idx) -= 2.0 * eps;
        const MatXd minus = eph.Eval(t_s, p);
        const double err_plus
            = (plus.leftCols(3) - baseline.leftCols(3)).rowwise().norm().maxCoeff();
        const double err_minus
            = (minus.leftCols(3) - baseline.leftCols(3)).rowwise().norm().maxCoeff();
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
             || name.rfind("M", 0) == 0 || name.rfind("u", 0) == 0;
    }
    bool IsEccentricityParam(const std::string& name) { return name.rfind("e", 0) == 0; }
    bool IsSemiMajorAxisParam(const std::string& name) { return name.rfind("a", 0) == 0; }

    template <typename EphemT>
    EphemerisWindowResult SummarizeWindows(const EphemT& eph, const std::vector<VecXd>& t_windows,
                                           const std::vector<MatXd>& rv_windows,
                                           double fit_window_min, double precision_m) {
      const int num_params = eph.NumParams();
      const std::vector<std::string> names = eph.ParamNames();
      const int nw = static_cast<int>(t_windows.size());
      VecXd max_k = VecXd::Zero(num_params);
      VecXd max_abs = VecXd::Zero(num_params);

      // Pass 1: fit every window and size each parameter's broadcast resolution
      // (the per-parameter fractional-bit count `max_k`, maximized across windows).
      std::vector<VecXd> params_per_window(nw);
      for (size_t w = 0; w < t_windows.size(); ++w) {
        params_per_window[w] = eph.Fit(t_windows[w], rv_windows[w]);
        for (int j = 2; j < num_params; ++j) {
          max_abs(j) = std::max(max_abs(j), std::abs(params_per_window[w](j)));
          const int k = SearchParamBits(eph, t_windows[w], params_per_window[w], j, precision_m);
          max_k(j) = std::max(max_k(j), static_cast<double>(k));
        }
      }

      // Pass 2: quantize each window's parameters to the chosen broadcast
      // resolution (LSB = 2^-max_k(j)) and evaluate the *realized* fit error, so
      // the reported accuracy is what a receiver decoding the quantized message
      // actually gets -- consistent with the `total_bits` datasize below. (t_ref /
      // t_fit, indices 0-1, are not broadcast-quantized.)
      double pos_rms_sum = 0.0, vel_rms_sum = 0.0, pos_p95_sum = 0.0, vel_p95_sum = 0.0;
      for (size_t w = 0; w < t_windows.size(); ++w) {
        VecXd qp = params_per_window[w];
        for (int j = 2; j < num_params; ++j) {
          const double lsb = std::pow(2.0, -max_k(j));
          qp(j) = std::round(qp(j) / lsb) * lsb;
        }
        const EphemerisFitErrorStats stats = eph.EvalError(t_windows[w], rv_windows[w], qp);
        pos_rms_sum += stats.rms_pos_m(3);
        vel_rms_sum += stats.rms_vel_mps(3);
        pos_p95_sum += stats.p95_pos_m(3);
        vel_p95_sum += stats.p95_vel_mps(3);
      }

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

    Frame ParseFrame(const Config& node, Frame fallback) {
      if (!node) return fallback;
      return enum_cast<Frame>(node.as<std::string>()).value();
    }

  }  // namespace

  // --- EphemerisSimulationConfig <- Config (YAML) translation -------------------

  EphemerisSimulationConfig ConfigToEphemerisConfig(Config& config) {
    EphemerisSimulationConfig c;
    c.start_epoch_utc = config["start_epoch_utc"].as<std::string>(c.start_epoch_utc);

    if (config["orbit"]) {
      Config o(config["orbit"]);
      c.orbit.a_m = o["a_m"].as<double>(c.orbit.a_m);
      c.orbit.ecc = o["ecc"].as<double>(c.orbit.ecc);
      c.orbit.inc_rad = o["inc_rad"].as<double>(c.orbit.inc_rad);
      c.orbit.raan_rad = o["raan_rad"].as<double>(c.orbit.raan_rad);
      c.orbit.argp_rad = o["argp_rad"].as<double>(c.orbit.argp_rad);
      c.orbit.m0_rad = o["m0_rad"].as<double>(c.orbit.m0_rad);
      c.orbit.coe_frame = ParseFrame(o["coe_frame"], c.orbit.coe_frame);
    }

    c.propagate_frame = ParseFrame(config["propagate_frame"], c.propagate_frame);
    c.duration_days = config["duration_days"].as<double>(c.duration_days);
    c.sample_dt_s = config["sample_dt_s"].as<double>(c.sample_dt_s);

    // Unified `force_model:` block (bodies list); falls back to the legacy scalar keys.
    if (config["force_model"]) {
      const ForceModelSpec fm = ParseForceModelSpec(config["force_model"]);
      c.moon_gravity_degree = fm.moon_degree;
      c.moon_gravity_order = fm.moon_order;
      c.include_earth = fm.include_earth;
      c.include_sun = fm.include_sun;
      c.use_relativity = fm.relativity;
    } else {
      c.moon_gravity_degree = config["moon_gravity_degree"].as<int>(c.moon_gravity_degree);
      c.moon_gravity_order = config["moon_gravity_order"].as<int>(c.moon_gravity_order);
      c.include_earth = config["include_earth"].as<bool>(c.include_earth);
      c.include_sun = config["include_sun"].as<bool>(c.include_sun);
      c.use_relativity = config["use_relativity"].as<bool>(c.use_relativity);
    }
    c.integration_step_s = config["integration_step_s"].as<double>(c.integration_step_s);

    if (config["fit_window_minutes"]) {
      c.fit_window_minutes.clear();
      for (const auto& item : config["fit_window_minutes"])
        c.fit_window_minutes.push_back(item.as<double>());
    }
    c.num_windows = config["num_windows"].as<int>(c.num_windows);

    c.cartesian_poly_order = config["cartesian_poly_order"].as<int>(c.cartesian_poly_order);
    c.cartesian_use_keplerian_baseline
        = config["cartesian_use_keplerian_baseline"].as<bool>(c.cartesian_use_keplerian_baseline);
    c.cartesian_num_fourier_terms
        = config["cartesian_num_fourier_terms"].as<int>(c.cartesian_num_fourier_terms);
    c.almanac_poly_order = config["almanac_poly_order"].as<int>(c.almanac_poly_order);
    c.almanac_num_fourier_terms
        = config["almanac_num_fourier_terms"].as<int>(c.almanac_num_fourier_terms);

    c.output_frame = ParseFrame(config["output_frame"], c.output_frame);
    c.datasize_precision_m = config["datasize_precision_m"].as<double>(c.datasize_precision_m);
    return c;
  }

  // --- EphemerisApp -------------------------------------------------------------

  EphemerisApp::EphemerisApp(Config& config) : Application(config) {
    cfg_ = ConfigToEphemerisConfig(config);
  }

  EphemerisApp::EphemerisApp(EphemerisSimulationConfig config) : cfg_(std::move(config)) {}

  void EphemerisApp::Setup() {
    LUPNT_CHECK(agent_, "Agent not set", "EphemerisApp");
    // The study is one-shot: schedule a single Step at t = 0. The heavy sweep runs
    // lazily on that first Step (so a programmatically-set config is honored).
    Simulation* sim = agent_->GetSimulation();
    sim->Schedule(
        0.0, [this](Real t) { Step(t); }, Event::SINGLE_EVENT, Event::Priority::APPLICATION);
  }

  void EphemerisApp::Step(Real /*t*/) {
    if (initialized_) return;
    Initialize();
    initialized_ = true;
  }

  void EphemerisApp::Initialize() {
    const EphemerisSimulationConfig& cfg = cfg_;
    LUPNT_CHECK(cfg.duration_days > 0.0, "duration_days must be positive", "EphemerisApp");
    LUPNT_CHECK(cfg.sample_dt_s > 0.0, "sample_dt_s must be positive", "EphemerisApp");
    LUPNT_CHECK(!cfg.fit_window_minutes.empty(), "fit_window_minutes must not be empty",
                "EphemerisApp");
    LUPNT_CHECK(cfg.num_windows >= 2, "num_windows must be at least 2", "EphemerisApp");
    LUPNT_CHECK(cfg.datasize_precision_m > 0.0, "datasize_precision_m must be positive",
                "EphemerisApp");

    // The truth trajectory is propagated in ABSOLUTE TDB seconds past J2000 (t0_tdb +
    // elapsed), passed directly to the orbit dynamics -- exactly as the former monolithic
    // driver did. Orbit propagation adds `GetLupntEpoch()` internally, so we reset the
    // global epoch to 0 here (the base `Simulation` sets it from the scenario `epoch:`) to
    // avoid double-counting the epoch and to reproduce the monolith's numerics bit-for-bit.
    SetLupntEpoch(0.0);

    const Real t0_tdb = ConvertTime(GregorianToTime(cfg.start_epoch_utc), Time::UTC, Time::TDB);

    Vec6 coe0;
    coe0 << cfg.orbit.a_m, cfg.orbit.ecc, cfg.orbit.inc_rad, cfg.orbit.raan_rad, cfg.orbit.argp_rad,
        cfg.orbit.m0_rad;
    const Vec6 rv0_coe_frame = ClassicalToCart(coe0, Real(GM_MOON));
    const Vec6 rv0
        = cfg.orbit.coe_frame == cfg.propagate_frame
              ? rv0_coe_frame
              : ConvertFrame(t0_tdb, rv0_coe_frame, cfg.orbit.coe_frame, cfg.propagate_frame);

    Ptr<NBodyDynamics> dynamics = MakeTruthDynamics(cfg);

    const int n
        = static_cast<int>(std::floor(cfg.duration_days * SECS_DAY / cfg.sample_dt_s + 1.0e-9)) + 1;
    VecX tspan(n);
    for (int i = 0; i < n; ++i) tspan(i) = i * cfg.sample_dt_s;
    VecX tfs(n);
    for (int i = 0; i < n; ++i) tfs(i) = t0_tdb + tspan(i);

    Cart6 x0(rv0, cfg.propagate_frame);
    MatX rv_prop = dynamics->Propagate(x0, tfs);

    // Fit/broadcast in the configured output frame: convert the truth trajectory
    // from the (inertial) propagate frame to output_frame per epoch if they differ
    // (e.g. MOON_CI -> MOON_PA for a rotating Moon-fixed broadcast frame).
    if (cfg.output_frame != cfg.propagate_frame) {
      for (int i = 0; i < n; ++i) {
        rv_prop.row(i) = ConvertFrame(tfs(i), Vec6(rv_prop.row(i).transpose()), cfg.propagate_frame,
                                      cfg.output_frame)
                             .transpose();
      }
    }

    results_.t_truth_s = tspan.cast<double>();
    results_.rv_truth = rv_prop.cast<double>();
    results_.cartesian_results.clear();
    results_.almanac_results.clear();

    EphemerisFitOptions eph_opts;
    eph_opts.poly_order = cfg.cartesian_poly_order;
    eph_opts.use_keplerian_baseline = cfg.cartesian_use_keplerian_baseline;
    eph_opts.gm = GM_MOON;
    eph_opts.frame = cfg.output_frame;
    eph_opts.num_fourier_terms = cfg.cartesian_num_fourier_terms;
    const CartesianEphemeris cart_eph(eph_opts);

    AlmanacFitOptions alm_opts;
    alm_opts.poly_order = cfg.almanac_poly_order;
    alm_opts.num_fourier_terms = cfg.almanac_num_fourier_terms;
    alm_opts.gm = GM_MOON;
    alm_opts.frame = cfg.output_frame;
    const Almanac almanac(alm_opts);

    auto pbar = Logger::GetProgressBar(static_cast<int>(cfg.fit_window_minutes.size()),
                                       "Ephemeris datasize/accuracy sweep");
    int prog = 0;
    for (double fit_min : cfg.fit_window_minutes) {
      const double window_s = fit_min * SECS_MINUTE;
      const int window_points = static_cast<int>(std::round(window_s / cfg.sample_dt_s)) + 1;
      LUPNT_CHECK(window_points >= 2 && window_points <= n,
                  "fit_window_minutes must fit within the sampled truth trajectory",
                  "EphemerisApp");

      const int margin = window_points / 2;
      const int span = std::max(1, n - 1 - 2 * margin);
      std::vector<VecXd> t_windows;
      std::vector<MatXd> rv_windows;
      for (int w = 0; w < cfg.num_windows; ++w) {
        const double frac
            = cfg.num_windows > 1 ? static_cast<double>(w) / (cfg.num_windows - 1) : 0.5;
        int start = margin + static_cast<int>(std::round(frac * span)) - window_points / 2;
        start = std::max(0, std::min(start, n - window_points));
        t_windows.push_back(results_.t_truth_s.segment(start, window_points));
        rv_windows.push_back(results_.rv_truth.middleRows(start, window_points));
      }

      results_.cartesian_results.push_back(
          SummarizeWindows(cart_eph, t_windows, rv_windows, fit_min, cfg.datasize_precision_m));
      results_.almanac_results.push_back(
          SummarizeWindows(almanac, t_windows, rv_windows, fit_min, cfg.datasize_precision_m));

      pbar->Update(++prog);
    }
    pbar->Finish();
  }

  REGISTER_FACTORY_CLASS(Application, EphemerisApp)

}  // namespace lupnt
