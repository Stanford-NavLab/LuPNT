#include "lupnt/measurements/gnss_measurement.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "lupnt/attitude/gnss_attitude.h"
#include "lupnt/attitude/gnss_yaw_steering.h"
#include "lupnt/conversions/frame_converter.h"
#include "lupnt/conversions/time_conversions.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/error.h"
#include "lupnt/devices/gnss_device.h"
#include "lupnt/measurements/comms_utils.h"
#include "lupnt/numerics/interpolation.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {

  namespace {
    RowVec3d ToRowVec3d(const Vec3& v) {
      RowVec3d row;
      row << v(0).val(), v(1).val(), v(2).val();
      return row;
    }

    Vec3d ToVec3d(const Vec3& v) { return Vec3d(v(0).val(), v(1).val(), v(2).val()); }

    Real ObservableSigma(const GnssChannel& channel, GnssObservable observable) {
      switch (observable) {
        case GnssObservable::PSEUDORANGE: return channel.sigma_pseudorange_m;
        case GnssObservable::DOPPLER: return channel.sigma_doppler_hz;
        case GnssObservable::CARRIER_PHASE: return channel.sigma_carrier_phase_cycles;
        default: LUPNT_CHECK(false, "Unknown GNSS observable", "GnssMeasurement");
      }
    }

    Real OneIfFinitePositive(Real value) {
      return std::isfinite(value.val()) && value.val() > 0.0 ? value : Real(1.0);
    }

    Real ConvertGnssEpoch(Real t, Time from, Time to) {
      if (from == to) return t;
      if (IsCoordinateTimeScale(from) && IsCoordinateTimeScale(to)) {
        return ConvertCoordinateTime(t, from, to);
      }
      return ConvertTime(t, from, to);
    }

    Real TransmitterRelativisticCorrectionSeconds(const Vec6& rv_tx) {
      return -2.0 * rv_tx.head(3).dot(rv_tx.tail(3)) / (Real(C) * Real(C));
    }

    Vec6 ConvertEciTransmitStateToMeasurementFrame(const Vec6& rv_tx_eci, Real t_ephem,
                                                   const GnssMeasurementOptions& options) {
      if (options.frame == Frame::ECI) return rv_tx_eci;
      Real t_tdb = ConvertGnssEpoch(t_ephem, options.ephemeris_time_scale, Time::TDB);
      return ConvertFrame(t_tdb, rv_tx_eci, Frame::ECI, options.frame);
    }

    Vec6 ResolveTransmitState(const GnssChannel& channel, const Vec3& r_rx,
                              const GnssMeasurementOptions& options) {
      if (!options.recompute_transmit_time_from_ephemeris || !channel.HasEphemeris()) {
        return channel.GetTransmitState(channel.transmit_time);
      }

      Epoch t_rx_ephem = channel.receive_time.To(channel.ephemeris_time_scale);
      // Solve for the light-time OFFSET dtau = t_tx - t_rx (a small, ~0.08 s
      // full-precision quantity), holding the absolute receive epoch fixed.
      // Iterating on t_tx = t_rx - rho/C directly would subtract rho/C from an
      // ~8e8 s epoch and lose ~1e-7 s (~30 m) of light time to the double ULP;
      // carrying the offset avoids that cancellation entirely. The absolute
      // transmit epoch (t_rx + dtau) is only ever used to *evaluate* the
      // ephemeris, where the residual ULP maps to sub-mm position error.
      Real dtau = Real(0.0);
      Vec6 rv_tx = channel.GetTransmitState(t_rx_ephem);
      if (!options.solve_light_time) return rv_tx;

      for (int iter = 0; iter < options.light_time_max_iterations; iter++) {
        Real rho = (r_rx - rv_tx.head(3)).norm();
        Real dtau_next = -rho / Real(C);
        bool converged = abs(dtau_next - dtau) < options.light_time_tolerance_s;
        dtau = dtau_next;
        rv_tx = channel.GetTransmitState(t_rx_ephem + dtau);
        if (converged) break;
      }
      return rv_tx;
    }

    Real LinkBudget(Real P_tx_dbw, Real G_tx_db, Real G_rx_db, Real range_m, GnssFreq freq,
                    const GnssReceiverParams& rx_params) {
      constexpr double k_B = 1.38e-23;  // [J/K] Boltzmann constant

      auto it = GNSS_FREQ_MAP.find(freq);
      LUPNT_CHECK(it != GNSS_FREQ_MAP.end(), "Unknown GNSS frequency", "GNSSMeasurements");
      Real freq_Hz = it->second;

      Real L_fs = FreeSpacePathLoss(range_m, freq_Hz);
      double kb_db = 10.0 * std::log10(k_B);

      return P_tx_dbw + G_tx_db + G_rx_db - rx_params.L_atm - L_fs - rx_params.L_ad
             - rx_params.L_pol - kb_db - 10.0 * log10(rx_params.T_eff);
    }

    double GnssFrequencyHz(GnssFreq freq) {
      auto it = GNSS_FREQ_MAP.find(freq);
      LUPNT_CHECK(it != GNSS_FREQ_MAP.end(), "Unknown GNSS frequency", "GNSSMeasurements");
      return it->second.val();
    }
  }  // namespace

  Real GnssChannel::Wavelength() const {
    auto it = GNSS_FREQ_MAP.find(frequency);
    LUPNT_CHECK(it != GNSS_FREQ_MAP.end(), "Unknown GNSS frequency", "GnssChannel");
    return Real(C) / it->second;
  }

  bool GnssChannel::HasEphemeris() const {
    if (!ephemeris_chebyshev.segments.empty()) return true;
    return ephemeris_times.size() > 0 && ephemeris_tx_states.rows() == ephemeris_times.size()
           && ephemeris_tx_states.cols() == 6;
  }

  Vec6 GnssChannel::GetTransmitState(const Epoch& t) const {
    // Cached-state fast path: return the stored state when `t` is the same
    // instant as `transmit_time`.
    //
    // With `Epoch` this comparison is *exact* -- the split representation
    // resolves far below the ~0.25 us ULP of an absolute seconds-from-J2000
    // double. The previous `Real` version needed an epoch-magnitude-scaled
    // tolerance, because any fixed tolerance below one ULP could only ever fire
    // on bit-equality, making the branch rounding-dependent. Carrying the scale
    // in the value also means a same-numeric-value-different-scale epoch no
    // longer aliases onto the cache.
    if (t == transmit_time) return tx_state;
    if (!HasEphemeris()) return tx_state;

    // The ephemeris tables are numeric and indexed in `ephemeris_time_scale`,
    // so collapse to seconds here. That collapse is lossy (it reintroduces the
    // epoch ULP), but only for the *evaluation argument*: at ~km/s transmitter
    // speeds a 0.25 us epoch error maps to sub-millimetre position error. The
    // light time itself is carried as an offset by the caller and never goes
    // through this step.
    const Real t_sec = t.ToSeconds();

    if (!ephemeris_chebyshev.segments.empty()) {
      VecX state_x;
      LUPNT_CHECK(ephemeris_chebyshev.Eval(t_sec, &state_x, nullptr),
                  "GNSS channel transmit time is outside the Chebyshev ephemeris window",
                  "GnssChannel");
      return Vec6(state_x);
    }

    Vec6 state;
    for (int k = 0; k < 6; k++) {
      VecXd col = ephemeris_tx_states.col(k);
      state(k) = LinearInterp1d(ephemeris_times, col, t_sec.val());
    }
    return state;
  }

  Real GnssChannel::EffectiveTransmitClockBiasSeconds() const {
    return tx_clock_bias_s + relativistic_correction_s - group_delay_s;
  }

  Real GnssChannel::EffectiveTransmitClockDrift() const { return tx_clock_drift; }

  VecXd GnssMeasurementValue::AsVector(const std::vector<GnssObservable>& observables) const {
    VecXd y(observables.size());
    for (int i = 0; i < static_cast<int>(observables.size()); i++) {
      switch (observables[i]) {
        case GnssObservable::PSEUDORANGE: y(i) = pseudorange_m.val(); break;
        case GnssObservable::DOPPLER: y(i) = doppler_hz.val(); break;
        case GnssObservable::CARRIER_PHASE:
          y(i) = (carrier_phase_cycles + carrier_integer_cycles).val();
          break;
        default: LUPNT_CHECK(false, "Unknown GNSS observable", "GnssMeasurementValue");
      }
    }
    return y;
  }

  GnssMeasurement::GnssMeasurement(const GnssChannel& channel) : channel_(channel) {}

  GnssMeasurement::GnssMeasurement(const GnssChannel& channel, GnssMeasurementOptions options)
      : channel_(channel), options_(std::move(options)) {}

  Real GnssMeasurement::ClockBiasToMeters(Real bias, ClockBiasUnit unit) {
    return ClockDynamics::BiasUnitsToSeconds(bias, unit) * C;
  }

  Real GnssMeasurement::ClockDriftToMetersPerSecond(Real drift, ClockBiasUnit unit) {
    return ClockDynamics::BiasUnitsToSeconds(drift, unit) * C;
  }

  Real GnssMeasurement::ClockBiasMetersDerivative(ClockBiasUnit unit) {
    return Real(C / ClockDynamics::SecondsToBiasUnitScale(unit));
  }

  Real GnssMeasurement::ClockDriftMetersPerSecondDerivative(ClockBiasUnit unit) {
    return Real(C / ClockDynamics::SecondsToBiasUnitScale(unit));
  }

  GnssMeasurementValue GnssMeasurement::ComputeValue(const State& user_state,
                                                     const GnssMeasurementOptions& options) const {
    const auto& idx = options.indices;
    LUPNT_CHECK(user_state.size() >= idx.position + 3, "User state does not contain position",
                "GnssMeasurement");
    LUPNT_CHECK(user_state.size() >= idx.velocity + 3, "User state does not contain velocity",
                "GnssMeasurement");

    Vec3 r_rx = user_state.segment(idx.position, 3);
    Vec3 v_rx = user_state.segment(idx.velocity, 3);
    Vec6 tx_state = ResolveTransmitState(channel_, r_rx, options);
    Vec3 r_tx = tx_state.head(3);
    Vec3 v_tx = tx_state.tail(3);

    Vec3 dr = r_rx - r_tx;
    Vec3 dv = v_rx - v_tx;
    Real rho = dr.norm();
    LUPNT_CHECK(rho > 0.0, "GNSS range is zero", "GnssMeasurement");
    Real rho_dot = dr.dot(dv) / rho;

    Real rx_bias_m = 0.0;
    Real rx_drift_mps = 0.0;
    if (idx.clock_bias >= 0 && idx.clock_bias < user_state.size()) {
      rx_bias_m = ClockBiasToMeters(user_state(idx.clock_bias), options.clock_bias_unit);
    }
    if (idx.clock_drift >= 0 && idx.clock_drift < user_state.size()) {
      rx_drift_mps
          = ClockDriftToMetersPerSecond(user_state(idx.clock_drift), options.clock_bias_unit);
    }

    Real tx_bias_m = C * channel_.EffectiveTransmitClockBiasSeconds();
    Real tx_drift_mps = C * channel_.EffectiveTransmitClockDrift();
    Real code_delay_m = channel_.shapiro_delay_m + channel_.ionosphere_plasma_delay_m;
    Real carrier_delay_m = channel_.shapiro_delay_m - channel_.ionosphere_plasma_delay_m;
    Real range_with_clock = rho + code_delay_m + rx_bias_m - tx_bias_m;
    Real range_rate_with_clock = rho_dot + rx_drift_mps - tx_drift_mps;
    Real lambda = channel_.Wavelength();

    Real integer_cycles = channel_.integer_ambiguity_cycles;
    if (idx.carrier_integer >= 0 && idx.carrier_integer < user_state.size()) {
      integer_cycles = user_state(idx.carrier_integer);
    }

    GnssMeasurementValue value;
    value.pseudorange_m = range_with_clock;
    value.doppler_hz = -range_rate_with_clock / lambda;
    value.carrier_phase_cycles
        = (rho + carrier_delay_m + rx_bias_m - tx_bias_m) / lambda + channel_.phase_bias_cycles;
    value.carrier_integer_cycles = integer_cycles;
    return value;
  }

  VecXd GnssMeasurement::ComputeVector(const State& user_state, MatXd* H,
                                       const GnssMeasurementOptions& options) const {
    GnssMeasurementValue value = ComputeValue(user_state, options);
    VecXd y = value.AsVector(options.observables);

    if (H == nullptr) return y;

    const auto& idx = options.indices;
    H->setZero(options.observables.size(), user_state.size());

    Vec3 r_rx = user_state.segment(idx.position, 3);
    Vec3 v_rx = user_state.segment(idx.velocity, 3);
    Vec6 tx_state = ResolveTransmitState(channel_, r_rx, options);
    Vec3 r_tx = tx_state.head(3);
    Vec3 v_tx = tx_state.tail(3);
    Vec3 dr = r_rx - r_tx;
    Vec3 dv = v_rx - v_tx;
    Real rho = dr.norm();
    Real rho_dot = dr.dot(dv) / rho;
    Vec3 u = dr / rho;
    Vec3 d_rhodot_dr = dv / rho - rho_dot * dr / (rho * rho);
    Vec3 d_rhodot_dv = u;
    Real lambda = channel_.Wavelength();
    Real db_m = ClockBiasMetersDerivative(options.clock_bias_unit);
    Real dd_mps = ClockDriftMetersPerSecondDerivative(options.clock_bias_unit);

    for (int row = 0; row < static_cast<int>(options.observables.size()); row++) {
      switch (options.observables[row]) {
        case GnssObservable::PSEUDORANGE:
          H->block(row, idx.position, 1, 3) = ToRowVec3d(u);
          if (idx.clock_bias >= 0 && idx.clock_bias < user_state.size()) {
            (*H)(row, idx.clock_bias) = db_m.val();
          }
          break;

        case GnssObservable::DOPPLER:
          H->block(row, idx.position, 1, 3) = ToRowVec3d(-d_rhodot_dr / lambda);
          H->block(row, idx.velocity, 1, 3) = ToRowVec3d(-d_rhodot_dv / lambda);
          if (idx.clock_drift >= 0 && idx.clock_drift < user_state.size()) {
            (*H)(row, idx.clock_drift) = (-dd_mps / lambda).val();
          }
          break;

        case GnssObservable::CARRIER_PHASE:
          H->block(row, idx.position, 1, 3) = ToRowVec3d(u / lambda);
          if (idx.clock_bias >= 0 && idx.clock_bias < user_state.size()) {
            (*H)(row, idx.clock_bias) = (db_m / lambda).val();
          }
          if (idx.carrier_integer >= 0 && idx.carrier_integer < user_state.size()) {
            (*H)(row, idx.carrier_integer) = 1.0;
          }
          break;

        default: LUPNT_CHECK(false, "Unknown GNSS observable", "GnssMeasurement");
      }
    }

    return y;
  }

  MeasData GnssMeasurement::Compute(const State& x, MatXd* H) const {
    MeasData md;
    md.timestamp = channel_.receive_time.ToSeconds();
    md.value = ComputeVector(x, H, options_);
    md.covariance = MatXd::Zero(md.value.size(), md.value.size());
    for (int i = 0; i < md.value.size(); i++) {
      Real sigma = OneIfFinitePositive(ObservableSigma(channel_, options_.observables[i]));
      md.covariance(i, i) = (sigma * sigma).val();
    }
    return md;
  }

  std::vector<GNSSMeasurementsEpoch> GnssMeasurement::Precompute(
      const GNSSMeasurements& measurements, const std::vector<Real>& receive_times,
      const std::vector<State>& user_states, bool compute_jacobians) {
    LUPNT_CHECK(receive_times.size() == user_states.size(),
                "Receive times and user states must have same size", "GnssMeasurement");
    std::vector<GNSSMeasurementsEpoch> epochs;
    epochs.reserve(receive_times.size());
    for (size_t i = 0; i < receive_times.size(); i++) {
      MatXd H;
      epochs.push_back(
          measurements.Compute(receive_times[i], user_states[i], compute_jacobians ? &H : nullptr));
    }
    return epochs;
  }

  GNSSMeasurements::GNSSMeasurements(Ptr<GnssConstellation> constellation) {
    constellations_.emplace_back(constellation, frequency_);
  }

  void GNSSMeasurements::AddConstellation(Ptr<GnssConstellation> constellation,
                                          GnssFreq frequency) {
    constellations_.emplace_back(constellation, frequency);
  }

  Vec3 GNSSMeasurements::SunPosition(Real t) const {
    if (sun_position_provider_) return sun_position_provider_(t);
    return Vec3(0.0, AU, 0.0);
  }

  Vec3 GNSSMeasurements::BoresightTarget(Real t) const {
    if (boresight_target_provider_) return boresight_target_provider_(t);
    (void)t;
    return Vec3::Zero();
  }

  GnssChannel GNSSMeasurements::BuildChannelForPrn(int prn, Real receive_time,
                                                   const State& user_state,
                                                   bool compute_ionosphere_plasma_delay,
                                                   const Ptr<GnssConstellation>& constellation,
                                                   GnssFreq frequency) const {
    LUPNT_CHECK(options_.ephemeris_time_scale == Time::TAI,
                "GnssConstellation ephemerides are stored and interpolated in TAI seconds",
                "GNSSMeasurements");
    LUPNT_CHECK(options_.frame != Frame::UNDEFINED,
                "GNSSMeasurements requires receiver and transmitter states in a defined common "
                "inertial frame",
                "GNSSMeasurements");

    const auto& idx = options_.indices;
    Vec3 r_rx = user_state.segment(idx.position, 3);

    Epoch t_rx_ephem = Epoch::FromSeconds(receive_time, options_.receive_time_scale)
                           .To(options_.ephemeris_time_scale);
    // Solve for the light-time OFFSET dtau = t_tx - t_rx (small, full precision)
    // with the absolute receive epoch fixed. Carried on an `Epoch` the offset
    // lives in the fractional second, so it is never absorbed into the epoch's
    // ~0.25 us ULP the way `t_rx - rho/C` on a bare double would be (~30 m of
    // light time). The absolute transmit epoch is formed only to evaluate the
    // ephemeris, where the residual maps to sub-mm position error.
    Real dtau = Real(0.0);
    Epoch t_tx_ephem = t_rx_ephem;
    Vec6 rv_tx = ConvertEciTransmitStateToMeasurementFrame(
        constellation->GetSatelliteStateEci(prn, t_tx_ephem.ToSeconds()), t_tx_ephem.ToSeconds(),
        options_);

    if (options_.solve_light_time) {
      for (int iter = 0; iter < options_.light_time_max_iterations; iter++) {
        Real rho = (r_rx - rv_tx.head(3)).norm();
        Real dtau_next = -rho / Real(C);
        bool converged = abs(dtau_next - dtau) < options_.light_time_tolerance_s;
        dtau = dtau_next;
        t_tx_ephem = t_rx_ephem + dtau;
        rv_tx = ConvertEciTransmitStateToMeasurementFrame(
            constellation->GetSatelliteStateEci(prn, t_tx_ephem.ToSeconds()),
            t_tx_ephem.ToSeconds(), options_);
        if (converged) break;
      }
    }

    GnssChannel channel;
    channel.gnss_const = constellation->GetGnssConst();
    channel.prn = prn;
    channel.frequency = frequency;
    channel.receive_time = Epoch::FromSeconds(receive_time, options_.receive_time_scale);
    channel.transmit_time = t_tx_ephem;
    channel.ephemeris_time_scale = options_.ephemeris_time_scale;
    channel.frame = options_.frame;
    channel.tx_state = rv_tx;
    channel.ephemeris_times = constellation->GetEphemerisTimesTai();
    if (options_.frame == Frame::ECI) {
      channel.ephemeris_tx_states = constellation->GetSatelliteStateHistoryEci(prn);
      channel.ephemeris_chebyshev = constellation->GetSatelliteStateChebyshevEci(prn);
    } else {
      const VecXd& ephemeris_times = channel.ephemeris_times;
      const MatXd& rv_eci_history = constellation->GetSatelliteStateHistoryEci(prn);
      channel.ephemeris_tx_states.resize(rv_eci_history.rows(), rv_eci_history.cols());
      for (int i = 0; i < rv_eci_history.rows(); ++i) {
        Vec6 rv_eci = rv_eci_history.row(i).transpose();
        Vec6 rv_frame
            = ConvertEciTransmitStateToMeasurementFrame(rv_eci, ephemeris_times(i), options_);
        channel.ephemeris_tx_states.row(i) = rv_frame.transpose();
      }
    }
    if (options_.apply_transmitter_relativity) {
      channel.relativistic_correction_s = TransmitterRelativisticCorrectionSeconds(rv_tx);
    }
    channel.shapiro_delay_m = ComputeShapiroDelay(receive_time, r_rx, rv_tx.head(3));
    if (compute_ionosphere_plasma_delay) {
      channel.ionosphere_plasma_delay_m
          = ComputeIonospherePlasmaDelay(receive_time, r_rx, rv_tx.head(3), frequency);
    }
    return channel;
  }

  Real GNSSMeasurements::ComputeShapiroDelay(Real receive_time, const Vec3& r_rx_eci,
                                             const Vec3& r_tx_eci) const {
    if (!options_.apply_shapiro_delay) return Real(0.0);

    Vec3 r_sun = SunPosition(receive_time);
    Vec3 r_rx = r_rx_eci - r_sun;
    Vec3 r_tx = r_tx_eci - r_sun;
    Real rx_norm = r_rx.norm();
    Real tx_norm = r_tx.norm();
    Real rho = (r_rx_eci - r_tx_eci).norm();
    if (rx_norm <= 0.0 || tx_norm <= 0.0) return Real(0.0);

    Real numerator = rx_norm + tx_norm + rho;
    Real denominator = rx_norm + tx_norm - rho;
    if (denominator <= 0.0 || numerator <= 0.0) return Real(0.0);

    return 2.0 * Real(GM_SUN) / (Real(C) * Real(C)) * log(numerator / denominator);
  }

  Real GNSSMeasurements::ComputeIonospherePlasmaDelay(Real receive_time, const Vec3& r_rx_eci,
                                                      const Vec3& r_tx_eci, GnssFreq freq) const {
    if (!options_.apply_ionosphere_plasma_delay) return Real(0.0);
    if (custom_ionosphere_plasma_delay_model_) {
      return custom_ionosphere_plasma_delay_model_(receive_time, r_rx_eci, r_tx_eci, freq);
    }
    if (ionosphere_plasma_raytrace_options_) {
      return ComputeIonospherePlasmaRayTraceDelay(receive_time, r_rx_eci, r_tx_eci, freq);
    }
    return options_.default_ionosphere_plasma_delay_m;
  }

  Real GNSSMeasurements::ComputeIonospherePlasmaRayTraceDelay(Real receive_time, const Vec3& r_rx,
                                                              const Vec3& r_tx,
                                                              GnssFreq freq) const {
    LUPNT_CHECK(ionosphere_plasma_raytrace_options_.has_value(),
                "Ionosphere/plasma raytrace options are not set", "GNSSMeasurements");
    const auto& ray_opts = *ionosphere_plasma_raytrace_options_;

    pecsim::RayTraceConfig config = ray_opts.config;
    if (ray_opts.use_channel_frequency) config.freq_Hz = GnssFrequencyHz(freq);
    LUPNT_CHECK(std::isfinite(config.freq_Hz) && config.freq_Hz > 0.0,
                "Raytrace frequency must be positive", "GNSSMeasurements");

    if (ray_opts.iri_model != pecsim::IRIModel::NONE) {
      pecsim::set_iri_model(ray_opts.iri_model);
    }

    Real t_raytrace = ConvertGnssEpoch(receive_time, options_.receive_time_scale,
                                       ray_opts.raytrace_epoch_scale);
    Real t_tdb = ConvertGnssEpoch(receive_time, options_.receive_time_scale, Time::TDB);

    Vec3 r_rx_ray = r_rx;
    Vec3 r_tx_ray = r_tx;
    if (ray_opts.raytrace_frame != options_.frame) {
      r_rx_ray = ConvertFrame(t_tdb, r_rx, options_.frame, ray_opts.raytrace_frame);
      r_tx_ray = ConvertFrame(t_tdb, r_tx, options_.frame, ray_opts.raytrace_frame);
    }

    pecsim::PathProfile profile = pecsim::trace_ray(
        t_raytrace.val(), ToVec3d(r_tx_ray) / 1000.0, ToVec3d(r_rx_ray) / 1000.0, config,
        ray_opts.debug_propagation, ray_opts.debug_correction);

    Real delay_m = profile.tec_delay_m;
    if (ray_opts.delay_mode == GnssIonospherePlasmaRayTraceDelayMode::TEC_PLUS_HIGHER_ORDER) {
      delay_m += profile.second_delay_m + profile.third_delay_m;
    }
    return delay_m;
  }

  Real GNSSMeasurements::ComputeCN0(const GnssChannel& channel, const Vec3& r_rx_eci,
                                    Real receive_time,
                                    const Ptr<GnssConstellation>& constellation) const {
    if (!constellation->HasTransmitterInfo(channel.prn, channel.frequency)) return Real(NAN);

    Vec6 rv_tx = channel.GetTransmitState(channel.transmit_time);
    Vec3 r_tx = rv_tx.head(3);
    Vec3 dr = r_rx_eci - r_tx;
    Real range = dr.norm();
    Vec3 u_tx2rx = dr.normalized();
    Vec3 u_rx2tx = -u_tx2rx;

    Vec3 r_boresight_target_eci = BoresightTarget(receive_time);
    Vec3 u_rx2body = (r_boresight_target_eci - r_rx_eci).normalized();

    Real t_tdb = ConvertGnssEpoch(receive_time, options_.receive_time_scale, Time::TDB);
    Vec6 rv_tx_gcrf = ConvertFrame(t_tdb, rv_tx, options_.frame, Frame::GCRF);
    Vec3 r_rx_gcrf = ConvertFrame(t_tdb, r_rx_eci, options_.frame, Frame::GCRF);
    Vec3 r_sun_gcrf = ConvertFrame(t_tdb, SunPosition(receive_time), options_.frame, Frame::GCRF);
    Vec3 u_tx2rx_gcrf = (r_rx_gcrf - rv_tx_gcrf.head(3)).normalized();

    Vec3 ex, ey, ez;
    if (options_.tx_yaw_model == GnssMeasurementOptions::TxYawModel::DEDICATED) {
      // Block-specific dedicated eclipse yaw law, evaluated continuously from geometry.
      const Vec3 r_tx_g = rv_tx_gcrf.head(3);
      const Vec3 v_tx_g = rv_tx_gcrf.tail(3);
      Real beta = GnssYawSteering::BetaAngle(r_tx_g, v_tx_g, r_sun_gcrf);
      Real mu = GnssYawSteering::OrbitAngle(r_tx_g, v_tx_g, r_sun_gcrf);
      Real phi_nom = GnssYawSteering::NominalYawAngle(beta, mu);
      Real phi_ded;
      switch (channel.gnss_const) {
        case GnssConst::GALILEO:
          phi_ded = GnssYawSteering::GalileoIovEclipseYawAngle(GnssYawSteering::OrbitNoonAngle(mu),
                                                               beta, phi_nom);
          break;
        case GnssConst::GPS:
          phi_ded = GnssYawSteering::Gps3EclipseYawAngle(beta, mu, phi_nom);
          break;
        default:
          phi_ded = phi_nom;  // no dedicated law for this system -> nominal
          break;
      }
      GnssAttitude::ComputeFromYawAngle(r_tx_g, v_tx_g, phi_ded, ex, ey, ez);
    } else {
      GnssAttitude::Compute(rv_tx_gcrf.head(3), Vec3(rv_tx_gcrf.tail(3)), r_sun_gcrf, ex, ey, ez);
    }
    Real phi_tx = safe_acos(u_tx2rx_gcrf.dot(ez));
    Real theta_tx = atan2(u_tx2rx_gcrf.dot(ey), u_tx2rx_gcrf.dot(ex));
    Real phi_rx = safe_acos(u_rx2body.dot(u_rx2tx));

    const Antenna& tx_antenna
        = constellation->GetTransmitterAntenna(channel.prn, channel.frequency);
    Real G_tx = (options_.tx_gain_model == GnssMeasurementOptions::TxGainModel::AZIMUTH_AVERAGED)
                    ? tx_antenna.ComputeGainAzimuthAveraged(phi_tx)
                    : tx_antenna.ComputeGain(theta_tx, phi_tx);
    Real G_rx = rx_antenna_.ComputeGain(0.0, phi_rx);
    Real P_tx = constellation->GetTransmitPowerDbw(channel.prn, channel.frequency);

    return LinkBudget(P_tx, G_tx, G_rx, range, channel.frequency, rx_params_);
  }

  Real GNSSMeasurements::ComputeSigmaRange(Real cn0_dbhz, GnssFreq freq) const {
    Real CN0_w = DecibelToDecimal(cn0_dbhz);
    auto it = GNSS_RC_MAP.find(freq);
    LUPNT_CHECK(it != GNSS_RC_MAP.end(), "Unknown GNSS frequency", "GNSSMeasurements");
    Real Rc = it->second;
    Real Tc = 1.0 / Rc;

    DllParams params;
    params.B_dll = rx_params_.Bn;
    params.T_i = rx_params_.T;
    params.d = rx_params_.D;
    params.Tc = Tc;
    params.B_fe = rx_params_.b * Rc;

    // DLL (code-tracking) thermal noise, root-sum-squared with the optional
    // broadcast ephemeris/clock pseudorange error terms (Mina et al. 2025, Eq. 50).
    Real sigma_dll = SigmaDll(params, CN0_w) * (Real(C) * Tc);
    Real sigma_eph = rx_params_.sigma_pr_eph_m;
    Real sigma_clk = rx_params_.sigma_pr_clk_m;
    return sqrt(sigma_dll * sigma_dll + sigma_eph * sigma_eph + sigma_clk * sigma_clk);
  }

  Real GNSSMeasurements::ComputeSigmaRangeRate(Real cn0_dbhz, GnssFreq freq) const {
    Real CN0_w = DecibelToDecimal(cn0_dbhz);
    auto it = GNSS_FREQ_MAP.find(freq);
    LUPNT_CHECK(it != GNSS_FREQ_MAP.end(), "Unknown GNSS frequency", "GNSSMeasurements");
    Real lambda = Real(C) / it->second;

    return lambda / (2.0 * PI * rx_params_.T)
           * sqrt(4.0 * rx_params_.Bf / CN0_w * (1.0 + 1.0 / (rx_params_.T * CN0_w)));
  }

  Real GNSSMeasurements::ComputeSigmaCarrierPhase(Real cn0_dbhz, GnssFreq freq) const {
    Real CN0_w = DecibelToDecimal(cn0_dbhz);
    auto it = GNSS_FREQ_MAP.find(freq);
    LUPNT_CHECK(it != GNSS_FREQ_MAP.end(), "Unknown GNSS frequency", "GNSSMeasurements");
    Real lambda = Real(C) / it->second;

    PllParams params;
    params.B_pll = rx_params_.Bp;
    params.T_i = rx_params_.T;

    // PLL (carrier-tracking) thermal noise, root-sum-squared with the optional
    // oscillator Allan-deviation and vibration-induced phase-noise terms
    // (Mina et al. 2025, Eqs. 54-55). All extra terms default to 0.
    Real sigma_pll_m = lambda / (2.0 * PI) * SigmaPll(params, CN0_w);
    // Second-order-PLL Allan-deviation phase error, expressed in meters (Eq. 55).
    Real sigma_ad_m = rx_params_.allan_deviation > 0.0
                          ? Real(C) / (360.0 * 144.0) * rx_params_.allan_deviation / rx_params_.Bp
                          : Real(0.0);
    Real sigma_vib_m = lambda / 360.0 * rx_params_.sigma_vib_deg;
    return sqrt(sigma_pll_m * sigma_pll_m + sigma_ad_m * sigma_ad_m + sigma_vib_m * sigma_vib_m);
  }

  std::vector<GnssChannel> GNSSMeasurements::BuildChannels(Real receive_time,
                                                           const State& user_state) const {
    return BuildChannels(receive_time, user_state, true);
  }

  std::vector<GnssChannel> GNSSMeasurements::BuildChannels(
      Real receive_time, const State& user_state, bool compute_ionosphere_plasma_delay) const {
    LUPNT_CHECK(!constellations_.empty(), "No GNSS constellations added", "GNSSMeasurements");
    LUPNT_CHECK(user_state.size() >= options_.indices.position + 3,
                "User state does not contain position", "GNSSMeasurements");
    LUPNT_CHECK(user_state.size() >= options_.indices.velocity + 3,
                "User state does not contain velocity", "GNSSMeasurements");

    std::vector<GnssChannel> channels;
    Vec3 r_rx = user_state.segment(options_.indices.position, 3);

    // Resolve per-epoch body positions once (avoid redundant callbacks per PRN).
    std::vector<Vec3> body_positions;
    body_positions.reserve(occluding_bodies_.size());
    for (const auto& body : occluding_bodies_) {
      body_positions.push_back(body.position_provider ? body.position_provider(receive_time)
                                                      : body.position_m);
    }

    for (const auto& [constellation, frequency] : constellations_) {
      LUPNT_CHECK(constellation != nullptr, "GNSS constellation is not set", "GNSSMeasurements");

      for (int prn : constellation->GetPrns()) {
        if (constellation->IsFaultPrn(prn)) continue;

        GnssChannel channel
            = BuildChannelForPrn(prn, receive_time, user_state, compute_ionosphere_plasma_delay,
                                 constellation, frequency);
        Vec3 r_tx = channel.tx_state.head(3);

        if (options_.apply_visibility) {
          bool visible = true;
          for (size_t bi = 0; bi < occluding_bodies_.size(); ++bi) {
            if (!ComputeVisibility(r_rx, r_tx, occluding_bodies_[bi].radius_m,
                                   body_positions[bi])) {
              visible = false;
              break;
            }
          }
          if (!visible) continue;
        }

        channel.cn0_dbhz = ComputeCN0(channel, r_rx, receive_time, constellation);
        if (options_.apply_cn0_threshold) {
          bool in_lock = tracking_prns_.count({channel.gnss_const, prn}) > 0;
          Real threshold = in_lock ? options_.cn0_tracking_threshold_dbhz
                                   : options_.cn0_acquisition_threshold_dbhz;
          if (!(channel.cn0_dbhz > threshold)) continue;
        }

        if (std::isfinite(channel.cn0_dbhz.val())) {
          Real sigma_range = ComputeSigmaRange(channel.cn0_dbhz, channel.frequency);
          Real sigma_rate = ComputeSigmaRangeRate(channel.cn0_dbhz, channel.frequency);
          Real sigma_phase = ComputeSigmaCarrierPhase(channel.cn0_dbhz, channel.frequency);
          channel.sigma_pseudorange_m = sigma_range;
          channel.sigma_doppler_hz = sigma_rate / channel.Wavelength();
          channel.sigma_carrier_phase_cycles = sigma_phase / channel.Wavelength();
        }

        channels.push_back(channel);
      }
    }

    // Update tracking state: satellites that survived this epoch stay in lock.
    tracking_prns_.clear();
    for (const auto& ch : channels) tracking_prns_.insert({ch.gnss_const, ch.prn});

    return channels;
  }

  GNSSMeasurementsEpoch GNSSMeasurements::ComputeFromChannels(Real receive_time,
                                                              const State& user_state,
                                                              std::vector<GnssChannel> channels,
                                                              MatXd* H) const {
    const int n_obs = static_cast<int>(options_.observables.size());
    const int n_rows = n_obs * static_cast<int>(channels.size());

    GNSSMeasurementsEpoch epoch;
    epoch.receive_time = Epoch::FromSeconds(receive_time, options_.receive_time_scale);
    epoch.channels = channels;
    epoch.values.resize(n_rows);
    epoch.jacobian.setZero(n_rows, user_state.size());
    epoch.covariance.setZero(n_rows, n_rows);

    for (int i = 0; i < static_cast<int>(channels.size()); i++) {
      GnssMeasurement measurement(channels[i]);
      MatXd H_i;
      VecXd y_i = measurement.ComputeVector(user_state, H != nullptr ? &H_i : nullptr, options_);
      epoch.values.segment(i * n_obs, n_obs) = y_i;
      if (H != nullptr) epoch.jacobian.block(i * n_obs, 0, n_obs, user_state.size()) = H_i;

      for (int j = 0; j < n_obs; j++) {
        Real sigma = OneIfFinitePositive(ObservableSigma(channels[i], options_.observables[j]));
        epoch.covariance(i * n_obs + j, i * n_obs + j) = (sigma * sigma).val();
      }
    }

    if (H != nullptr) *H = epoch.jacobian;
    return epoch;
  }

  GNSSMeasurementsEpoch GNSSMeasurements::Compute(Real receive_time, const State& user_state,
                                                  MatXd* H) const {
    return ComputeFromChannels(receive_time, user_state, BuildChannels(receive_time, user_state),
                               H);
  }

  std::vector<GNSSMeasurementsEpoch> GNSSMeasurements::Precompute(
      const std::vector<Real>& receive_times, const std::vector<State>& user_states,
      bool compute_jacobians) {
    LUPNT_CHECK(receive_times.size() == user_states.size(),
                "Receive times and user states must have same size", "GNSSMeasurements");
    precomputed_.clear();
    precomputed_.reserve(receive_times.size());

    if (options_.apply_ionosphere_plasma_delay && batch_custom_ionosphere_plasma_delay_model_) {
      std::vector<std::vector<GnssChannel>> channels_by_epoch;
      channels_by_epoch.reserve(receive_times.size());
      for (size_t i = 0; i < receive_times.size(); i++) {
        channels_by_epoch.push_back(BuildChannels(receive_times[i], user_states[i], false));
      }

      std::vector<std::vector<Real>> delays_m = batch_custom_ionosphere_plasma_delay_model_(
          receive_times, user_states, channels_by_epoch);
      LUPNT_CHECK(delays_m.size() == channels_by_epoch.size(),
                  "Batch ionosphere/plasma delay provider returned wrong number of epochs",
                  "GNSSMeasurements");

      for (size_t i = 0; i < channels_by_epoch.size(); i++) {
        LUPNT_CHECK(delays_m[i].size() == channels_by_epoch[i].size(),
                    "Batch ionosphere/plasma delay provider returned wrong number of channels",
                    "GNSSMeasurements");
        for (size_t j = 0; j < channels_by_epoch[i].size(); j++) {
          channels_by_epoch[i][j].ionosphere_plasma_delay_m = delays_m[i][j];
        }

        MatXd H;
        precomputed_.push_back(ComputeFromChannels(receive_times[i], user_states[i],
                                                   channels_by_epoch[i],
                                                   compute_jacobians ? &H : nullptr));
      }
      return precomputed_;
    }

    for (size_t i = 0; i < receive_times.size(); i++) {
      MatXd H;
      precomputed_.push_back(
          Compute(receive_times[i], user_states[i], compute_jacobians ? &H : nullptr));
    }
    return precomputed_;
  }

  FilterMeasurementFunction GNSSMeasurements::CreateFunction(Real t) const {
    return [manager = *this, t](const State& x, MatXd* H, MatXd* R) {
      GNSSMeasurementsEpoch epoch = manager.Compute(t, x, H);
      if (R != nullptr) *R = epoch.covariance;
      return epoch.values;
    };
  }

}  // namespace lupnt
