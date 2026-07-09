#pragma once

#include <functional>
#include <optional>
#include <set>
#include <utility>
#include <vector>

#include "lupnt/agents/gnss_constellation.h"
#include "lupnt/dynamics/clock_dynamics.h"
#include "lupnt/environment/plasma/gcpm/iri_interface.h"
#include "lupnt/environment/plasma/tec/raytrace.h"
#include "lupnt/measurements/measurement.h"
#include "lupnt/measurements/measurement_utils.h"
#include "lupnt/numerics/cheby_fit.h"
#include "lupnt/numerics/filters/filter.h"

namespace lupnt {

  enum class GnssObservable { PSEUDORANGE, DOPPLER, CARRIER_PHASE };

  struct GnssMeasurementStateIndices {
    int position = 0;
    int velocity = 3;
    int clock_bias = 6;
    int clock_drift = 7;
    int carrier_integer = -1;
  };

  struct GnssMeasurementOptions {
    std::vector<GnssObservable> observables
        = {GnssObservable::PSEUDORANGE, GnssObservable::DOPPLER, GnssObservable::CARRIER_PHASE};
    GnssMeasurementStateIndices indices;
    ClockBiasUnit clock_bias_unit = ClockBiasUnit::SECONDS;

    /// Epochs passed to `GNSSMeasurements::Compute` / `Precompute` are
    /// receiver signal-reception epochs in this time scale.
    Time receive_time_scale = Time::TAI;

    /// GNSS constellation ephemeris epochs are represented in this scale.
    /// `GnssConstellation` stores its tables as TAI seconds by convention.
    Time ephemeris_time_scale = Time::TAI;

    /// Receiver and transmitter states are expressed in this frame.
    Frame frame = Frame::ECI;

    bool solve_light_time = true;
    /// If true, `GnssMeasurement::Compute` resolves the transmit epoch from
    /// the current user state and the channel ephemeris coefficients instead
    /// of reusing the transmit epoch stored in the precomputed channel. This is
    /// intended for estimated/filter measurements; truth generation can keep
    /// the precomputed light-time solution tied to the truth state.
    bool recompute_transmit_time_from_ephemeris = false;
    int light_time_max_iterations = 10;
    Real light_time_tolerance_s = 1.0e-11;
    bool apply_transmitter_relativity = true;
    bool apply_shapiro_delay = true;
    /// Deprecated: Shapiro delay is currently evaluated for the Sun only.
    Real shapiro_mu = GM_SUN;
    /// Deprecated: Shapiro delay is currently evaluated using SunPosition().
    Vec3 shapiro_body_position = Vec3::Zero();
    bool apply_visibility = true;
    bool apply_cn0_threshold = true;
    Real cn0_threshold_dbhz = 15.0;              // Deprecated: use acq/tracking below
    Real cn0_acquisition_threshold_dbhz = 22.0;  // Min CN0 to acquire a new satellite [dBHz]
    Real cn0_tracking_threshold_dbhz = 20.0;     // Min CN0 to maintain an existing lock [dBHz]
    bool apply_ionosphere_plasma_delay = false;
    Real default_ionosphere_plasma_delay_m = 0.0;
  };

  struct GnssChannel {
    GnssConst gnss_const = GnssConst::GPS;
    int prn = 0;
    GnssFreq frequency = GnssFreq::L1;
    Real receive_time = 0.0;
    Time receive_time_scale = Time::TAI;
    Real transmit_time = 0.0;
    Time transmit_time_scale = Time::TAI;
    Time ephemeris_time_scale = Time::TAI;
    Frame frame = Frame::ECI;

    Vec6 tx_state = Vec6::Zero();
    VecXd ephemeris_times;
    MatXd ephemeris_tx_states;
    ChebyshevFitModel ephemeris_chebyshev;
    Real tx_clock_bias_s = 0.0;
    Real tx_clock_drift = 0.0;
    Real relativistic_correction_s = 0.0;
    Real group_delay_s = 0.0;
    Real shapiro_delay_m = 0.0;
    Real ionosphere_plasma_delay_m = 0.0;
    Real phase_bias_cycles = 0.0;
    Real integer_ambiguity_cycles = 0.0;

    Real cn0_dbhz = NAN;
    Real sigma_pseudorange_m = NAN;
    Real sigma_doppler_hz = NAN;
    Real sigma_carrier_phase_cycles = NAN;

    /// @brief Carrier wavelength `c / f` for this channel's GNSS frequency.
    ///
    /// Used by GnssMeasurement::ComputeValue/Compute to convert between range
    /// and Doppler/carrier-phase (cycles), and by GNSSMeasurements to convert
    /// the FLL/PLL range-rate and phase noise sigmas into
    /// `sigma_doppler_hz`/`sigma_carrier_phase_cycles`.
    ///
    /// @return Wavelength [m]; throws if `frequency` is not in `GNSS_FREQ_MAP`
    Real Wavelength() const;

    /// @brief True if this channel carries a transmitter ephemeris (either a
    /// fitted Chebyshev model or a sampled state history) usable by
    /// GetTransmitState(), as opposed to only the single `tx_state` snapshot.
    bool HasEphemeris() const;

    /// @brief Evaluate the transmitter (e.g. GNSS satellite) position/velocity
    /// state at epoch `t`.
    ///
    /// Used by GnssMeasurement::ComputeValue/Compute (via ResolveTransmitState)
    /// to get the transmitter state at the resolved transmit epoch when
    /// re-solving light time from the precomputed channel ephemeris. Falls
    /// back to the stored `tx_state` snapshot if HasEphemeris() is false.
    ///
    /// @param t Evaluation epoch, in `ephemeris_time_scale` [s]
    /// @return  Transmitter Cartesian state `[r; v]` [m, m/s] in `frame`
    Vec6 GetTransmitState(Real t) const;

    /// @brief Total transmitter clock bias to apply to the pseudorange,
    /// combining the raw clock bias, the special-relativistic periodic
    /// correction, and the broadcast group delay: `tx_clock_bias_s +
    /// relativistic_correction_s - group_delay_s`.
    ///
    /// Used by GnssMeasurement::ComputeValue to compute `tx_bias_m = C *
    /// EffectiveTransmitClockBiasSeconds()`, the transmitter-side clock term
    /// subtracted from the pseudorange/carrier-phase observables.
    ///
    /// @return Effective transmitter clock bias [s]
    Real EffectiveTransmitClockBiasSeconds() const;

    /// @brief Transmitter clock drift used for the Doppler/range-rate clock
    /// term (currently just `tx_clock_drift`).
    ///
    /// @return Effective transmitter clock drift [s/s]
    Real EffectiveTransmitClockDrift() const;
  };

  struct GnssMeasurementValue {
    Real pseudorange_m = NAN;
    Real doppler_hz = NAN;
    Real carrier_phase_cycles = NAN;
    Real carrier_integer_cycles = 0.0;

    /// @brief Stack the requested subset of this measurement's observables
    /// (pseudorange, Doppler, carrier phase) into a single vector, in the
    /// given order.
    ///
    /// Called by GnssMeasurement::Compute to assemble the measurement vector
    /// `y` returned to the EKF/UKF, matching the row ordering used when
    /// building the Jacobian `H`. For `CARRIER_PHASE`, the integer ambiguity
    /// (`carrier_integer_cycles`) is added to `carrier_phase_cycles`.
    ///
    /// @param observables Which observables to include, and in what order
    /// @return            Vector of observable values, one per entry of
    ///                     `observables` [m, Hz, or cycles depending on entry]
    VecXd AsVector(const std::vector<GnssObservable>& observables) const;
  };

  struct GNSSMeasurementsEpoch {
    /// Receiver signal-reception epoch.
    Real receive_time = 0.0;
    Time receive_time_scale = Time::TAI;
    Real time = 0.0;  // Backward-compatible alias of `receive_time`.

    std::vector<GnssChannel> channels;
    VecXd values;
    MatXd jacobian;
    MatXd covariance;
  };

  enum class GnssIonospherePlasmaRayTraceDelayMode { TEC_ONLY, TEC_PLUS_HIGHER_ORDER };

  struct GnssIonospherePlasmaRayTraceOptions {
    /// Ray-tracing configuration supplied by the caller. GNSSMeasurements does
    /// not fill in model parameters such as Kp, integrator, correction mode, or
    /// step sizes.
    pecsim::RayTraceConfig config;

    /// Frame expected by the plasma ray tracer. The raytrace implementation
    /// uses geocentric Earth-fixed coordinates by convention.
    Frame raytrace_frame = Frame::ECEF;

    /// Time scale passed to the plasma model. The raytrace routines expect
    /// seconds from J2000 in UTC by convention.
    Time raytrace_epoch_scale = Time::UTC;

    /// If true, overwrite `config.freq_Hz` with the GNSS channel frequency.
    bool use_channel_frequency = true;

    /// Optional IRI model selection. If left as NONE, the current global plasma
    /// IRI setting is used.
    pecsim::IRIModel iri_model = pecsim::IRIModel::NONE;

    GnssIonospherePlasmaRayTraceDelayMode delay_mode
        = GnssIonospherePlasmaRayTraceDelayMode::TEC_ONLY;
    bool debug_propagation = false;
    bool debug_correction = false;
  };

  class GNSSMeasurements;

  /// @brief Single-channel GNSS measurement model (`Measurement` subclass).
  ///
  /// Its `Config` is `GnssMeasurementOptions` (observable selection, state-index layout,
  /// clock-bias unit, frame, and light-time/relativity/Shapiro/plasma correction settings),
  /// stored as a member and consumed by the `Measurement` interface methods
  /// (`Compute`/`CreateFunction`). The explicit-config overloads
  /// `ComputeValue`/`ComputeVector` remain available for callers (e.g. the Lunar-GNSS ODTS
  /// combined/TDCP model) that vary the config per invocation.
  class GnssMeasurement : public MeasurementClone<GnssMeasurement> {
  public:
    /// @brief The measurement model's configuration type (observable selection, state
    /// mapping, correction settings). Alias of `GnssMeasurementOptions`.
    using Config = GnssMeasurementOptions;

    GnssMeasurement() = default;

    /// @brief Construct a single-channel GNSS measurement from a (typically
    /// precomputed) GnssChannel.
    ///
    /// Used by GNSSMeasurements::ComputeFromChannels to wrap each transmitter
    /// channel built by `BuildChannels` before computing its
    /// pseudorange/Doppler/carrier-phase observables.
    explicit GnssMeasurement(const GnssChannel& channel);

    /// @brief Construct from a channel and an explicit measurement config.
    GnssMeasurement(const GnssChannel& channel, GnssMeasurementOptions options);

    /// @brief Replace the GNSS channel (transmitter ephemeris, clock, and
    /// delay terms) used by ComputeValue()/ComputeVector().
    void SetChannel(const GnssChannel& channel) { channel_ = channel; }

    /// @brief Get the GNSS channel currently associated with this measurement.
    const GnssChannel& GetChannel() const { return channel_; }

    /// @brief Replace the measurement config used by the `Measurement` interface methods.
    void SetOptions(const GnssMeasurementOptions& options) { options_ = options; }

    /// @brief Get the measurement config currently in use.
    const GnssMeasurementOptions& GetOptions() const { return options_; }

    /// @brief Compute the pseudorange, Doppler, and carrier-phase values
    /// implied by a user (receiver) state and this measurement's channel.
    ///
    /// This is the core GNSS observation model used throughout the module:
    /// GnssMeasurement::Compute calls this to build the measurement vector,
    /// and GNSSMeasurements::ComputeFromChannels calls it per-channel when
    /// assembling a `GNSSMeasurementsEpoch`. It forms the geometric
    /// range/range-rate between `user_state` and the (light-time- and
    /// relativity-corrected) transmitter state, then applies the receiver and
    /// effective-transmitter clock bias/drift and the channel's
    /// Shapiro-delay/ionosphere-plasma corrections (added for the code/
    /// pseudorange observable, subtracted for the carrier-phase observable,
    /// per the standard dispersive-vs-non-dispersive sign convention).
    ///
    /// @param user_state Receiver state vector; must contain position and
    ///                    velocity at `options.indices.position`/`velocity`
    ///                    [m, m/s] in `options.frame`, and optionally receiver
    ///                    clock bias/drift and carrier-integer-ambiguity
    ///                    states at the other `indices` entries
    /// @param options    Observable selection, state-index layout, clock-bias
    ///                    unit, and light-time/relativity/Shapiro/plasma
    ///                    correction settings
    /// @return GnssMeasurementValue with pseudorange [m], Doppler [Hz], and
    ///         carrier-phase [cycles] (plus the carrier integer ambiguity)
    GnssMeasurementValue ComputeValue(const State& user_state,
                                      const GnssMeasurementOptions& options = {}) const;

    /// @brief Compute the GNSS measurement vector (and, optionally, its
    /// Jacobian) for the observables selected in `options`.
    ///
    /// The explicit-config workhorse behind the `Measurement` interface: it
    /// evaluates `y = h(x)` and `H = dh/dx` for one GNSS channel. Internally
    /// calls ComputeValue() for `y`, and (if `H` is requested) differentiates
    /// the range/range-rate/carrier-phase model analytically with respect to
    /// receiver position, velocity, clock bias/drift, and carrier-integer
    /// state, per the observables in `options.observables`.
    ///
    /// @param user_state Receiver state vector, as in ComputeValue()
    /// @param[out] H      If non-null, filled with the
    ///                     `observables.size() x user_state.size()` Jacobian
    ///                     of `y` with respect to `user_state`
    /// @param options    Same options as ComputeValue()
    /// @return Measurement vector `y`, one entry per `options.observables`
    ///         [m, Hz, or cycles depending on entry]
    VecXd ComputeVector(const State& user_state, MatXd* H = nullptr,
                        const GnssMeasurementOptions& options = {}) const;

    /// @brief `Measurement` interface: predict `z = h(x)`, `R`, and (if `H != nullptr`) the
    /// analytic Jacobian `H = dh/dx`, using the stored config `GetOptions()`.
    ///
    /// The diagonal `R` covariance is built from `GetChannel()`'s
    /// `sigma_pseudorange_m`/`sigma_doppler_hz`/`sigma_carrier_phase_cycles` (each defaulted
    /// to 1 if not finite/positive). `CreateFunction()` (inherited) wraps this into a
    /// `FilterMeasurementFunction` for direct use with the `Filter` family.
    MeasData Compute(const State& x, MatXd* H = nullptr) const override;

    /// @brief Batch-compute a `GNSSMeasurementsEpoch` (channels + values [+
    /// Jacobians]) for each given receive time / user state pair.
    ///
    /// Thin wrapper around `GNSSMeasurements::Compute` used to generate a
    /// time series of GNSS measurement epochs (e.g. for truth-data generation
    /// or batch filter processing) without manually looping over
    /// `receive_times`/`user_states`.
    ///
    /// @param measurements    Configured GNSSMeasurements (constellation,
    ///                          options, antennas, etc.) used to build channels
    /// @param receive_times   Receiver signal-reception epochs
    ///                          [s, `measurements.GetOptions().receive_time_scale`]
    /// @param user_states     Receiver states at each `receive_times` entry,
    ///                          one-to-one with `receive_times`
    /// @param compute_jacobians If true, also compute and store each epoch's
    ///                          measurement Jacobian
    /// @return One `GNSSMeasurementsEpoch` per `receive_times` entry
    static std::vector<GNSSMeasurementsEpoch> Precompute(const GNSSMeasurements& measurements,
                                                         const std::vector<Real>& receive_times,
                                                         const std::vector<State>& user_states,
                                                         bool compute_jacobians = false);

  private:
    GnssChannel channel_;
    GnssMeasurementOptions options_;  ///< config consumed by the `Measurement` interface

    /// @brief Convert a receiver clock bias from `unit` to an equivalent
    /// pseudorange offset, via `ClockDynamics::BiasUnitsToSeconds(bias, unit) * C`.
    ///
    /// Used by ComputeValue() to fold the receiver clock-bias state
    /// (`user_state(idx.clock_bias)`, stored in `options.clock_bias_unit`)
    /// into the pseudorange/carrier-phase observables in meters.
    ///
    /// @param bias Receiver clock bias in `unit`
    /// @param unit Unit the clock-bias state is stored in (seconds, meters, or kilometers)
    /// @return     Equivalent range bias [m]
    static Real ClockBiasToMeters(Real bias, ClockBiasUnit unit);

    /// @brief Convert a receiver clock drift from `unit`/s to an equivalent
    /// range-rate offset, via `ClockDynamics::BiasUnitsToSeconds(drift, unit) * C`.
    ///
    /// Used by ComputeValue() to fold the receiver clock-drift state into the
    /// Doppler observable in meters per second.
    ///
    /// @param drift Receiver clock drift in `unit` per second
    /// @param unit  Unit the clock-drift state is stored in
    /// @return      Equivalent range-rate bias [m/s]
    static Real ClockDriftToMetersPerSecond(Real drift, ClockBiasUnit unit);

    /// @brief Partial derivative of ClockBiasToMeters() with respect to the
    /// receiver clock-bias state, i.e. `C / SecondsToBiasUnitScale(unit)`.
    ///
    /// Used by Compute() to populate the clock-bias column of the
    /// pseudorange/carrier-phase rows of the Jacobian `H`.
    ///
    /// @param unit Unit the clock-bias state is stored in
    /// @return     `d(range bias [m]) / d(clock-bias state [unit])`
    static Real ClockBiasMetersDerivative(ClockBiasUnit unit);

    /// @brief Partial derivative of ClockDriftToMetersPerSecond() with
    /// respect to the receiver clock-drift state (same scale factor as
    /// ClockBiasMetersDerivative()).
    ///
    /// Used by Compute() to populate the clock-drift column of the Doppler
    /// row of the Jacobian `H`.
    ///
    /// @param unit Unit the clock-drift state is stored in
    /// @return     `d(range-rate bias [m/s]) / d(clock-drift state [unit/s])`
    static Real ClockDriftMetersPerSecondDerivative(ClockBiasUnit unit);
  };

  class GNSSMeasurements {
  public:
    using VectorProvider = std::function<Vec3(Real)>;
    using CustomIonospherePlasmaDelayModel
        = std::function<Real(Real, const Vec3&, const Vec3&, GnssFreq)>;
    using BatchCustomIonospherePlasmaDelayModel = std::function<std::vector<std::vector<Real>>(
        const std::vector<Real>&, const std::vector<State>&,
        const std::vector<std::vector<GnssChannel>>&)>;

    GNSSMeasurements() = default;

    /// @brief Construct with a single GNSS constellation on L1 (convenience overload).
    explicit GNSSMeasurements(Ptr<GnssConstellation> constellation);

    /// @brief Add a constellation/frequency pair to process in BuildChannels().
    ///
    /// Channels from all added (constellation, frequency) pairs are concatenated.
    void AddConstellation(Ptr<GnssConstellation> constellation, GnssFreq frequency);

    /// @brief Replace the first (or only) constellation; clears all others.
    void SetConstellation(Ptr<GnssConstellation> constellation) {
      constellations_.clear();
      constellations_.emplace_back(constellation, frequency_);
    }

    /// @brief Get the first constellation (backward-compat accessor).
    Ptr<GnssConstellation> GetConstellation() const {
      return constellations_.empty() ? nullptr : constellations_[0].first;
    }

    /// @brief Set the frequency for all current constellations (and as default for future ones).
    void SetFrequency(GnssFreq frequency) {
      frequency_ = frequency;
      for (auto& p : constellations_) p.second = frequency;
    }

    /// @brief Get the frequency of the first constellation (backward-compat accessor).
    GnssFreq GetFrequency() const {
      return constellations_.empty() ? frequency_ : constellations_[0].second;
    }

    /// @brief Replace the measurement options (observable selection, state
    /// layout, light-time/relativity/Shapiro/plasma settings, etc.) used by
    /// BuildChannels()/Compute().
    void SetOptions(const GnssMeasurementOptions& options) { options_ = options; }

    /// @brief Get the measurement options currently in use.
    const GnssMeasurementOptions& GetOptions() const { return options_; }

    /// @brief Set the callback used by SunPosition() to get the Sun's
    /// position at a given epoch (used for Shapiro delay and GNSS transmitter
    /// attitude/antenna-gain computations).
    ///
    /// If unset, SunPosition() falls back to a fixed `(0, AU, 0)` placeholder.
    void SetSunPositionProvider(VectorProvider provider) { sun_position_provider_ = provider; }

    /// @brief Set the callback used by BoresightTarget() to get the
    /// receiver's nadir/boresight-pointing target position (e.g. the Moon or
    /// Earth center), used by ComputeCN0() to evaluate the receiver antenna's
    /// off-boresight angle to the transmitter.
    void SetBoresightTargetProvider(VectorProvider provider) {
      boresight_target_provider_ = provider;
    }

    /// @brief Set the spherical bodies (e.g. the Moon or Earth) checked for
    /// line-of-sight occlusion by BuildChannels() when
    /// `options.apply_visibility` is true (via ComputeVisibility()).
    void SetOccludingBodies(const std::vector<GnssOccludingBody>& bodies) {
      occluding_bodies_ = bodies;
    }

    /// @brief Set the receiver tracking-loop parameters (DLL/PLL/FLL
    /// bandwidths, integration time, etc.) used by ComputeSigmaRange(),
    /// ComputeSigmaRangeRate(), and ComputeSigmaCarrierPhase() to derive
    /// per-channel measurement noise sigmas from CN0.
    void SetReceiverParams(const GnssReceiverParams& params) { rx_params_ = params; }

    /// @brief Set the receiver antenna gain pattern used by ComputeCN0() to
    /// evaluate `G_rx` in the link budget.
    void SetReceiverAntenna(const Antenna& antenna) { rx_antenna_ = antenna; }

    /// @brief Set the CN0 threshold [dB-Hz] below which a channel is dropped
    /// by BuildChannels() when `options.apply_cn0_threshold` is true.
    /// Sets both acquisition and tracking thresholds to the same value.
    void SetCN0Threshold(Real cn0_threshold_dbhz) {
      options_.cn0_threshold_dbhz = cn0_threshold_dbhz;
      options_.cn0_acquisition_threshold_dbhz = cn0_threshold_dbhz;
      options_.cn0_tracking_threshold_dbhz = cn0_threshold_dbhz;
    }

    /// @brief Reset the internal tracking state, forcing all satellites to
    /// re-acquire at the acquisition threshold on the next BuildChannels() call.
    void ResetTracking() { tracking_prns_.clear(); }

    /// @brief Set a custom per-channel ionosphere/plasma delay model,
    /// overriding the ray-trace and default-constant delay options.
    ///
    /// If set, ComputeIonospherePlasmaDelay() calls this model instead of
    /// using `ionosphere_plasma_raytrace_options_` or
    /// `options_.default_ionosphere_plasma_delay_m`, taking the receive time
    /// and receiver/transmitter ECI positions and returning the (dispersive)
    /// delay for a given GNSS frequency.
    void SetCustomIonospherePlasmaDelayModel(CustomIonospherePlasmaDelayModel model) {
      custom_ionosphere_plasma_delay_model_ = model;
    }

    /// @brief Set a custom batch ionosphere/plasma delay model evaluated once
    /// per epoch for all channels, used by Precompute() instead of
    /// per-channel calls to ComputeIonospherePlasmaDelay() when
    /// `options.apply_ionosphere_plasma_delay` is true.
    void SetBatchCustomIonospherePlasmaDelayModel(BatchCustomIonospherePlasmaDelayModel model) {
      batch_custom_ionosphere_plasma_delay_model_ = model;
    }

    /// @brief Deprecated alias for SetCustomIonospherePlasmaDelayModel().
    [[deprecated("Use SetCustomIonospherePlasmaDelayModel")]]
    void SetIonospherePlasmaDelayProvider(CustomIonospherePlasmaDelayModel provider) {
      SetCustomIonospherePlasmaDelayModel(provider);
    }

    /// @brief Deprecated alias for SetBatchCustomIonospherePlasmaDelayModel().
    [[deprecated("Use SetBatchCustomIonospherePlasmaDelayModel")]]
    void SetBatchIonospherePlasmaDelayProvider(BatchCustomIonospherePlasmaDelayModel provider) {
      SetBatchCustomIonospherePlasmaDelayModel(provider);
    }

    /// @brief Configure ray-traced ionosphere/plasma delay computation (used
    /// by ComputeIonospherePlasmaRayTraceDelay() via
    /// ComputeIonospherePlasmaDelay()) when no custom model is set.
    void SetIonospherePlasmaRayTraceOptions(const GnssIonospherePlasmaRayTraceOptions& options) {
      ionosphere_plasma_raytrace_options_ = options;
    }

    /// @brief Disable ray-traced ionosphere/plasma delay computation,
    /// reverting ComputeIonospherePlasmaDelay() to the custom model (if any)
    /// or `options_.default_ionosphere_plasma_delay_m`.
    void ClearIonospherePlasmaRayTraceOptions() { ionosphere_plasma_raytrace_options_.reset(); }

    /// @brief Build the set of visible, above-CN0-threshold GNSS channels for
    /// every PRN in the constellation at a receiver epoch/state.
    ///
    /// This is the main entry point used by Compute()/Precompute() (and, via
    /// GnssMeasurement::Precompute, by truth/filter pipelines) to turn a
    /// receiver state into the list of GnssChannel transmitters that should
    /// contribute measurements at this epoch: for each non-faulted PRN it
    /// resolves the light-time-corrected transmitter state
    /// (BuildChannelForPrn), applies the visibility test (ComputeVisibility),
    /// computes CN0 (ComputeCN0) and drops channels below
    /// `options_.cn0_threshold_dbhz` if `options_.apply_cn0_threshold` is
    /// set, and fills in the per-channel measurement-noise sigmas.
    ///
    /// @param receive_time Receiver signal-reception epoch
    ///                       [s, `options_.receive_time_scale`]
    /// @param user_state   Receiver state; must contain position/velocity at
    ///                       `options_.indices.position`/`velocity`
    ///                       [m, m/s] in `options_.frame`
    /// @return One GnssChannel per visible, sufficiently-strong transmitter
    std::vector<GnssChannel> BuildChannels(Real receive_time, const State& user_state) const;

    /// @brief Build channels and compute the full GNSS measurement epoch
    /// (values, and optionally the Jacobian/covariance) for a receiver
    /// epoch/state.
    ///
    /// This is the primary per-epoch API used by the EKF/UKF measurement
    /// update (directly or via CreateFunction()) and by GnssMeasurement::Precompute:
    /// it calls BuildChannels() to get the visible transmitters, then
    /// ComputeFromChannels() to evaluate each channel's
    /// pseudorange/Doppler/carrier-phase observables.
    ///
    /// @param receive_time Receiver signal-reception epoch
    ///                       [s, `options_.receive_time_scale`]
    /// @param user_state   Receiver state, as in BuildChannels()
    /// @param[out] H        If non-null, filled with the stacked measurement
    ///                       Jacobian (same as `epoch.jacobian`)
    /// @return GNSSMeasurementsEpoch containing the built channels, stacked
    ///         measurement vector, Jacobian, and diagonal covariance
    GNSSMeasurementsEpoch Compute(Real receive_time, const State& user_state,
                                  MatXd* H = nullptr) const;

    /// @brief Compute and cache a GNSSMeasurementsEpoch for each of a series
    /// of receiver epochs/states.
    ///
    /// Used to generate a full simulation-length sequence of GNSS measurement
    /// epochs (e.g. for truth-data generation ahead of a filter run). If
    /// `options_.apply_ionosphere_plasma_delay` is set and a batch custom
    /// ionosphere/plasma delay model has been installed
    /// (SetBatchCustomIonospherePlasmaDelayModel), channels for all epochs are
    /// built first (without per-channel delays), the batch model is evaluated
    /// once for all epochs/channels, and the resulting delays are applied
    /// before computing each epoch's measurement values; otherwise each epoch
    /// is computed independently via Compute(). Results are also stored in
    /// `precomputed_` (see GetPrecomputed()).
    ///
    /// @param receive_times   Receiver signal-reception epochs
    ///                          [s, `options_.receive_time_scale`]
    /// @param user_states     Receiver states, one-to-one with `receive_times`
    /// @param compute_jacobians If true, also compute and store each epoch's
    ///                          measurement Jacobian
    /// @return One GNSSMeasurementsEpoch per `receive_times` entry
    std::vector<GNSSMeasurementsEpoch> Precompute(const std::vector<Real>& receive_times,
                                                  const std::vector<State>& user_states,
                                                  bool compute_jacobians = false);

    /// @brief Get the measurement epochs computed by the most recent Precompute() call.
    const std::vector<GNSSMeasurementsEpoch>& GetPrecomputed() const { return precomputed_; }

    /// @brief Build a `FilterMeasurementFunction` that calls Compute() at a
    /// fixed receive time `t`.
    ///
    /// The returned closure captures a copy of this GNSSMeasurements
    /// (constellation, options, antennas, providers, etc.) and is suitable
    /// for direct use as a filter's measurement model at epoch `t`, returning
    /// the stacked measurement vector and (if requested) Jacobian/covariance
    /// for all visible channels.
    ///
    /// @param t Receiver signal-reception epoch [s, `options_.receive_time_scale`]
    /// @return Function `(x, H, R) -> y` computing the full GNSS measurement
    ///         epoch at `t` for receiver state `x`
    FilterMeasurementFunction CreateFunction(Real t) const;

  private:
    std::vector<std::pair<Ptr<GnssConstellation>, GnssFreq>> constellations_;
    GnssFreq frequency_ = GnssFreq::L1;
    // Satellites currently in lock: (GnssConst, PRN). Mutable so BuildChannels
    // (const) can update tracking state each epoch.
    mutable std::set<std::pair<GnssConst, int>> tracking_prns_;
    GnssMeasurementOptions options_;
    std::vector<GnssOccludingBody> occluding_bodies_;
    Antenna rx_antenna_;
    GnssReceiverParams rx_params_;
    VectorProvider sun_position_provider_;
    VectorProvider boresight_target_provider_;
    CustomIonospherePlasmaDelayModel custom_ionosphere_plasma_delay_model_;
    BatchCustomIonospherePlasmaDelayModel batch_custom_ionosphere_plasma_delay_model_;
    std::optional<GnssIonospherePlasmaRayTraceOptions> ionosphere_plasma_raytrace_options_;
    std::vector<GNSSMeasurementsEpoch> precomputed_;

    /// @brief Get the Sun's position at epoch `t`, via
    /// `sun_position_provider_` if set, otherwise a fixed `(0, AU, 0)` placeholder.
    ///
    /// Used by ComputeShapiroDelay() (Sun-mass gravitational delay) and
    /// ComputeCN0() (transmitter attitude/yaw-steering frame).
    ///
    /// @param t Epoch [s, `options_.receive_time_scale`]
    /// @return  Sun position [m] in `options_.frame`
    Vec3 SunPosition(Real t) const;

    /// @brief Get the receiver's nadir/boresight-pointing target position at
    /// epoch `t`, via `boresight_target_provider_` if set, otherwise the origin.
    ///
    /// Used by ComputeCN0() to determine the receiver antenna's
    /// off-boresight angle to the transmitter.
    ///
    /// @param t Epoch [s, `options_.receive_time_scale`]
    /// @return  Boresight target position [m] in `options_.frame`
    Vec3 BoresightTarget(Real t) const;

    /// @brief Build the GnssChannel for a single PRN at a receiver epoch/state.
    ///
    /// Called once per PRN by the private BuildChannels() overload. Resolves
    /// the light-time-corrected transmit epoch and transmitter state from the
    /// constellation ephemeris, fills in the channel's time scales/frame and
    /// ephemeris data (including the Chebyshev fit, for later
    /// re-interpolation by GnssChannel::GetTransmitState), applies the
    /// transmitter relativistic clock correction
    /// (`options_.apply_transmitter_relativity`), and computes the Shapiro
    /// delay and (if requested) ionosphere/plasma delay for this
    /// transmitter-receiver pair.
    ///
    /// @param prn          GNSS satellite PRN within `constellation_`
    /// @param receive_time Receiver signal-reception epoch
    ///                       [s, `options_.receive_time_scale`]
    /// @param user_state   Receiver state; must contain position at
    ///                       `options_.indices.position` [m] in `options_.frame`
    /// @param compute_ionosphere_plasma_delay If true, also populate
    ///        `channel.ionosphere_plasma_delay_m` via ComputeIonospherePlasmaDelay()
    /// @return Fully populated GnssChannel for `prn` (visibility/CN0 not yet applied)
    GnssChannel BuildChannelForPrn(int prn, Real receive_time, const State& user_state,
                                   bool compute_ionosphere_plasma_delay,
                                   const Ptr<GnssConstellation>& constellation,
                                   GnssFreq frequency) const;

    /// @brief BuildChannels() implementation with explicit control over
    /// per-channel ionosphere/plasma delay computation.
    ///
    /// The public BuildChannels(receive_time, user_state) overload delegates
    /// here with `compute_ionosphere_plasma_delay = true`. Precompute() calls
    /// this directly with `false` when a batch ionosphere/plasma delay model
    /// will fill in the delays afterward.
    ///
    /// @param receive_time Receiver signal-reception epoch
    ///                       [s, `options_.receive_time_scale`]
    /// @param user_state   Receiver state, as in the public BuildChannels()
    /// @param compute_ionosphere_plasma_delay Forwarded to BuildChannelForPrn()
    /// @return One GnssChannel per visible, sufficiently-strong transmitter
    std::vector<GnssChannel> BuildChannels(Real receive_time, const State& user_state,
                                           bool compute_ionosphere_plasma_delay) const;

    /// @brief Evaluate the GNSS measurement model for a given (already-built)
    /// set of channels and assemble them into a GNSSMeasurementsEpoch.
    ///
    /// Called by Compute() (with freshly-built channels) and by Precompute()
    /// (with channels that may have had ionosphere/plasma delays patched in
    /// by a batch model). For each channel, constructs a GnssMeasurement and
    /// calls GnssMeasurement::Compute to fill the corresponding rows of
    /// `epoch.values`/`epoch.jacobian`, and fills `epoch.covariance`'s
    /// diagonal from the channels' observable-noise sigmas.
    ///
    /// @param receive_time Receiver signal-reception epoch
    ///                       [s, `options_.receive_time_scale`]
    /// @param user_state   Receiver state, as in Compute()
    /// @param channels     Channels to evaluate (e.g. from BuildChannels())
    /// @param[out] H        If non-null, set to `epoch.jacobian`
    /// @return GNSSMeasurementsEpoch with `channels`, stacked `values`,
    ///         `jacobian`, and diagonal `covariance` filled in
    GNSSMeasurementsEpoch ComputeFromChannels(Real receive_time, const State& user_state,
                                              std::vector<GnssChannel> channels,
                                              MatXd* H = nullptr) const;

    /// @brief Compute the Sun-mass relativistic (Shapiro) signal delay
    /// between a receiver and transmitter position.
    ///
    /// Called by BuildChannelForPrn() to fill `channel.shapiro_delay_m` when
    /// `options_.apply_shapiro_delay` is set; this delay is added to the
    /// pseudorange/code observable and subtracted from the carrier-phase
    /// observable in GnssMeasurement::ComputeValue (per the
    /// non-dispersive-delay sign convention). Uses the standard
    /// general-relativistic light-bending formula
    /// `2*GM_sun/c^2 * ln((r_rx+r_tx+rho)/(r_rx+r_tx-rho))` with positions
    /// measured relative to the Sun (via SunPosition()).
    ///
    /// @param receive_time Receiver signal-reception epoch
    ///                       [s, `options_.receive_time_scale`]
    /// @param r_rx_eci     Receiver position [m] in `options_.frame`
    /// @param r_tx_eci     Transmitter position [m] in `options_.frame`
    /// @return Shapiro delay [m] (0 if disabled or geometrically degenerate)
    Real ComputeShapiroDelay(Real receive_time, const Vec3& r_rx_eci, const Vec3& r_tx_eci) const;

    /// @brief Compute the ionosphere/plasma (dispersive) signal delay for one
    /// transmitter-receiver pair and frequency.
    ///
    /// Called by BuildChannelForPrn() to fill
    /// `channel.ionosphere_plasma_delay_m` when
    /// `options_.apply_ionosphere_plasma_delay` is set. Dispatches, in order
    /// of priority: a custom model (`custom_ionosphere_plasma_delay_model_`,
    /// see SetCustomIonospherePlasmaDelayModel), a configured ray-trace model
    /// (ComputeIonospherePlasmaRayTraceDelay(), see
    /// SetIonospherePlasmaRayTraceOptions), or finally a constant
    /// (`options_.default_ionosphere_plasma_delay_m`).
    ///
    /// @param receive_time Receiver signal-reception epoch
    ///                       [s, `options_.receive_time_scale`]
    /// @param r_rx_eci     Receiver position [m] in `options_.frame`
    /// @param r_tx_eci     Transmitter position [m] in `options_.frame`
    /// @param freq         GNSS signal frequency the delay applies to
    /// @return Ionosphere/plasma delay [m] (0 if disabled)
    Real ComputeIonospherePlasmaDelay(Real receive_time, const Vec3& r_rx_eci, const Vec3& r_tx_eci,
                                      GnssFreq freq) const;

    /// @brief Compute the ionosphere/plasma delay by ray-tracing through a
    /// plasma/TEC model (GCPM/IRI) between the transmitter and receiver.
    ///
    /// Called by ComputeIonospherePlasmaDelay() when ray-trace options have
    /// been configured (SetIonospherePlasmaRayTraceOptions) and no custom
    /// model is set. Converts the receive epoch and receiver/transmitter
    /// positions into the ray tracer's expected time scale and frame
    /// (`raytrace_epoch_scale`/`raytrace_frame`), optionally overrides the
    /// configured frequency with `freq`, runs `pecsim::trace_ray`, and
    /// returns either the TEC-only delay or TEC plus higher-order terms
    /// depending on `delay_mode`.
    ///
    /// @param receive_time Receiver signal-reception epoch
    ///                       [s, `options_.receive_time_scale`]
    /// @param r_rx         Receiver position [m] in `options_.frame`
    /// @param r_tx         Transmitter position [m] in `options_.frame`
    /// @param freq         GNSS signal frequency the delay applies to
    /// @return Ray-traced ionosphere/plasma delay [m]
    Real ComputeIonospherePlasmaRayTraceDelay(Real receive_time, const Vec3& r_rx, const Vec3& r_tx,
                                              GnssFreq freq) const;

    /// @brief Estimate the carrier-to-noise density ratio (CN0) for a GNSS channel.
    ///
    /// Called by the private BuildChannels() overload after visibility
    /// filtering, to decide (via `options_.cn0_threshold_dbhz`) whether a
    /// channel is trackable and to drive the per-channel measurement-noise
    /// sigmas (ComputeSigmaRange/ComputeSigmaRangeRate/ComputeSigmaCarrierPhase).
    /// Computes the transmitter and receiver antenna gains
    /// (`G_tx` from the constellation's transmitter antenna pattern evaluated
    /// at the off-boresight/azimuth angles given the GNSS attitude frame,
    /// `G_rx` from `rx_antenna_` evaluated at the angle to BoresightTarget())
    /// and the transmit power, then evaluates the link budget at the current
    /// range. Returns NaN if the constellation has no transmitter info for
    /// this PRN/frequency.
    ///
    /// @param channel      Channel whose `prn`/`frequency`/`transmit_time`
    ///                       identify the transmitter
    /// @param r_rx_eci     Receiver position [m] in `options_.frame`
    /// @param receive_time Receiver signal-reception epoch
    ///                       [s, `options_.receive_time_scale`]
    /// @return Estimated CN0 [dB-Hz], or NaN if unavailable
    Real ComputeCN0(const GnssChannel& channel, const Vec3& r_rx_eci, Real receive_time,
                    const Ptr<GnssConstellation>& constellation) const;

    /// @brief Convert CN0 to a pseudorange (code-tracking) noise standard deviation.
    ///
    /// Called by the private BuildChannels() overload to fill
    /// `channel.sigma_pseudorange_m` for channels with a finite CN0.
    /// Builds DllParams from `rx_params_` and the frequency's chipping rate
    /// (`GNSS_RC_MAP`), and scales SigmaDll()'s chip-fraction error by `C * Tc`.
    ///
    /// @param cn0_dbhz CN0 [dB-Hz]
    /// @param freq     GNSS frequency (selects the spreading-code chipping rate)
    /// @return Pseudorange noise standard deviation [m]
    Real ComputeSigmaRange(Real cn0_dbhz, GnssFreq freq) const;

    /// @brief Convert CN0 to a range-rate (FLL) noise standard deviation,
    /// later divided by wavelength to get `sigma_doppler_hz`.
    ///
    /// Called by the private BuildChannels() overload to derive
    /// `channel.sigma_doppler_hz = ComputeSigmaRangeRate(...) /
    /// channel.Wavelength()` for channels with a finite CN0.
    ///
    /// @param cn0_dbhz CN0 [dB-Hz]
    /// @param freq     GNSS frequency (selects the carrier wavelength)
    /// @return Range-rate noise standard deviation [m/s]
    Real ComputeSigmaRangeRate(Real cn0_dbhz, GnssFreq freq) const;

    /// @brief Convert CN0 to a carrier-phase (PLL) noise standard deviation
    /// in meters, later divided by wavelength to get
    /// `sigma_carrier_phase_cycles`.
    ///
    /// Called by the private BuildChannels() overload to derive
    /// `channel.sigma_carrier_phase_cycles = ComputeSigmaCarrierPhase(...) /
    /// channel.Wavelength()` for channels with a finite CN0. Builds
    /// PllParams from `rx_params_` and scales SigmaPll()'s phase error [rad]
    /// to meters via `lambda / (2*pi)`.
    ///
    /// @param cn0_dbhz CN0 [dB-Hz]
    /// @param freq     GNSS frequency (selects the carrier wavelength)
    /// @return Carrier-phase noise standard deviation [m]
    Real ComputeSigmaCarrierPhase(Real cn0_dbhz, GnssFreq freq) const;
  };

}  // namespace lupnt
