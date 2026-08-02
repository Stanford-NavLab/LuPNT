#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  GnssChannel MakeSimpleChannel() {
    GnssChannel channel;
    channel.gnss_const = GnssConst::GPS;
    channel.prn = 1;
    channel.frequency = GnssFreq::L1;
    channel.frame = Frame::ECI;
    channel.tx_state << 0.0, 0.0, 0.0, 0.0, 100.0, 0.0;
    channel.tx_clock_bias_s = 2.0e-6;
    channel.tx_clock_drift = 3.0e-9;
    channel.integer_ambiguity_cycles = 12.0;
    channel.phase_bias_cycles = 0.25;
    channel.sigma_pseudorange_m = 2.0;
    channel.sigma_doppler_hz = 0.5;
    channel.sigma_carrier_phase_cycles = 0.01;
    return channel;
  }

  JointOrbitClockState MakeUserState(ClockBiasUnit unit) {
    Cart6 orbit(Vec6::Zero(), Frame::ECI);
    orbit.r() = Vec3(3.0, 4.0, 0.0);
    orbit.v() = Vec3(0.6, 0.8, 0.0);

    ClockState2 clock;
    clock.b() = ClockDynamics::SecondsToBiasUnits(1.0e-6, unit);
    clock.d() = ClockDynamics::SecondsToBiasUnits(2.0e-9, unit);
    clock.SetUnits(ClockDynamics::GetStateUnits(2, unit));

    return JointOrbitClockState(orbit, clock);
  }

  MatXd CircularOrbitEci(double radius_m, double phase0_rad, const VecXd& t_tai) {
    double omega = std::sqrt(GM_EARTH / (radius_m * radius_m * radius_m));
    MatXd rv(t_tai.size(), 6);
    for (int k = 0; k < t_tai.size(); k++) {
      double th = phase0_rad + omega * t_tai(k);
      rv(k, 0) = radius_m * std::cos(th);
      rv(k, 1) = radius_m * std::sin(th);
      rv(k, 2) = 0.0;
      rv(k, 3) = -radius_m * omega * std::sin(th);
      rv(k, 4) = radius_m * omega * std::cos(th);
      rv(k, 5) = 0.0;
    }
    return rv;
  }
}  // namespace

TEST_CASE("measurements.gnss_measurement.single_channel") {
  GnssChannel channel = MakeSimpleChannel();
  GnssMeasurement measurement(channel);

  SECTION("computes pseudorange, Doppler, and carrier phase with seconds clock state") {
    JointOrbitClockState user_state = MakeUserState(ClockBiasUnit::SECONDS);
    GnssMeasurementOptions options;
    options.clock_bias_unit = ClockBiasUnit::SECONDS;

    GnssMeasurementValue value = measurement.ComputeValue(user_state, options);
    Real lambda = channel.Wavelength();

    double expected_range = 5.0 + C * (1.0e-6 - 2.0e-6);
    double expected_rate = ((3.0 * 0.6) + (4.0 * (0.8 - 100.0))) / 5.0 + C * (2.0e-9 - 3.0e-9);
    double expected_doppler = -expected_rate / lambda.val();
    double expected_phase = expected_range / lambda.val() + 0.25;

    REQUIRE_THAT(value.pseudorange_m.val(), WithinAbs(expected_range, 1.0e-8));
    REQUIRE_THAT(value.doppler_hz.val(), WithinAbs(expected_doppler, 1.0e-8));
    REQUIRE_THAT(value.carrier_phase_cycles.val(), WithinAbs(expected_phase, 1.0e-8));
    REQUIRE_THAT(value.carrier_integer_cycles.val(), WithinAbs(12.0, 1.0e-12));

    VecXd y = measurement.ComputeVector(user_state, nullptr, options);
    REQUIRE(y.size() == 3);
    REQUIRE_THAT(y(0), WithinAbs(expected_range, 1.0e-8));
    REQUIRE_THAT(y(1), WithinAbs(expected_doppler, 1.0e-8));
    REQUIRE_THAT(y(2), WithinAbs(expected_phase + 12.0, 1.0e-8));
  }

  SECTION("respects meter-scaled clock states") {
    JointOrbitClockState user_state = MakeUserState(ClockBiasUnit::METERS);
    GnssMeasurementOptions options;
    options.clock_bias_unit = ClockBiasUnit::METERS;

    GnssMeasurementValue value = measurement.ComputeValue(user_state, options);
    double expected_range = 5.0 + C * (1.0e-6 - 2.0e-6);
    REQUIRE_THAT(value.pseudorange_m.val(), WithinAbs(expected_range, 1.0e-8));
  }

  SECTION("adds non-dispersive Shapiro and dispersive plasma corrections") {
    GnssChannel corrected_channel = channel;
    corrected_channel.shapiro_delay_m = 0.08;
    corrected_channel.ionosphere_plasma_delay_m = 3.5;
    GnssMeasurement corrected(corrected_channel);

    JointOrbitClockState user_state = MakeUserState(ClockBiasUnit::SECONDS);
    GnssMeasurementValue value = corrected.ComputeValue(user_state);
    Real lambda = corrected_channel.Wavelength();

    double base_range = 5.0 + C * (1.0e-6 - 2.0e-6);
    double expected_pseudorange = base_range + 0.08 + 3.5;
    double expected_carrier = (base_range + 0.08 - 3.5) / lambda.val() + 0.25;

    REQUIRE_THAT(value.pseudorange_m.val(), WithinAbs(expected_pseudorange, 1.0e-8));
    REQUIRE_THAT(value.carrier_phase_cycles.val(), WithinAbs(expected_carrier, 1.0e-8));
  }

  SECTION("uses an estimated carrier integer ambiguity when configured") {
    State user_state(9);
    user_state.setZero();
    user_state.segment(0, 6) = MakeUserState(ClockBiasUnit::SECONDS).head(6);
    user_state(6) = 1.0e-6;
    user_state(7) = 2.0e-9;
    user_state(8) = 42.0;

    GnssMeasurementOptions options;
    options.indices.carrier_integer = 8;

    GnssMeasurementValue value = measurement.ComputeValue(user_state, options);
    REQUIRE_THAT(value.carrier_integer_cycles.val(), WithinAbs(42.0, 1.0e-12));

    MatXd H;
    VecXd y = measurement.ComputeVector(user_state, &H, options);
    REQUIRE(y.size() == 3);
    REQUIRE(H.rows() == 3);
    REQUIRE(H.cols() == 9);
    REQUIRE_THAT(H(2, 8), WithinAbs(1.0, 1.0e-12));
  }

  SECTION("analytic Jacobian has expected geometry and clock scaling") {
    JointOrbitClockState user_state = MakeUserState(ClockBiasUnit::SECONDS);
    GnssMeasurementOptions options;
    options.clock_bias_unit = ClockBiasUnit::SECONDS;

    MatXd H;
    measurement.ComputeVector(user_state, &H, options);
    Real lambda = channel.Wavelength();

    REQUIRE(H.rows() == 3);
    REQUIRE(H.cols() == user_state.size());
    REQUIRE_THAT(H(0, 0), WithinAbs(0.6, 1.0e-12));
    REQUIRE_THAT(H(0, 1), WithinAbs(0.8, 1.0e-12));
    REQUIRE_THAT(H(0, 6), WithinAbs(C, 1.0e-6));
    REQUIRE_THAT(H(1, 3), WithinAbs(-0.6 / lambda.val(), 1.0e-12));
    REQUIRE_THAT(H(1, 4), WithinAbs(-0.8 / lambda.val(), 1.0e-12));
    REQUIRE_THAT(H(1, 7), WithinAbs(-C / lambda.val(), 1.0e-6));
    REQUIRE_THAT(H(2, 0), WithinAbs(0.6 / lambda.val(), 1.0e-12));
    REQUIRE_THAT(H(2, 1), WithinAbs(0.8 / lambda.val(), 1.0e-12));
  }

  SECTION("filter function supplies diagonal measurement covariance") {
    JointOrbitClockState user_state = MakeUserState(ClockBiasUnit::SECONDS);
    FilterMeasurementFunction f = measurement.CreateFunction();

    MatXd H;
    MatXd R;
    VecXd y = f(user_state, &H, &R);

    REQUIRE(y.size() == 3);
    REQUIRE(R.rows() == 3);
    REQUIRE_THAT(R(0, 0), WithinAbs(4.0, 1.0e-12));
    REQUIRE_THAT(R(1, 1), WithinAbs(0.25, 1.0e-12));
    REQUIRE_THAT(R(2, 2), WithinAbs(1.0e-4, 1.0e-12));
  }
}

TEST_CASE("measurements.gnss_measurements.manager") {
  VecXd t_ephem(2);
  t_ephem << -100.0, 3600.0;

  auto constellation = MakePtr<GnssConstellation>(GnssConst::GPS);
  constellation->SetSatelliteStates({1}, t_ephem, {CircularOrbitEci(26560e3, 0.0, t_ephem)});
  constellation->SetupTransmitters();

  Cart6 orbit(Vec6::Zero(), Frame::ECI);
  orbit.r() = Vec3(R_EARTH + 500e3, 0.0, 0.0);
  orbit.v() = Vec3::Zero();
  ClockState2 clock;
  JointOrbitClockState user_state(orbit, clock);

  GNSSMeasurements manager(constellation);
  manager.SetFrequency(GnssFreq::L1);
  GnssMeasurementOptions manager_options;
  manager_options.apply_cn0_threshold = false;
  manager.SetOptions(manager_options);
  manager.SetSunPositionProvider([](Real) { return Vec3(0.0, AU, 0.0); });
  manager.SetBoresightTargetProvider([](Real) { return Vec3::Zero(); });

  GnssOccludingBody earth;
  earth.radius_m = R_EARTH;
  earth.position_m = Vec3::Zero();
  manager.SetOccludingBodies({earth});

  MatXd H;
  GNSSMeasurementsEpoch epoch = manager.Compute(0.0, user_state, &H);

  REQUIRE(epoch.channels.size() == 1);
  REQUIRE(epoch.channels[0].prn == 1);
  REQUIRE(epoch.receive_time.scale() == Time::TAI);
  REQUIRE(epoch.channels[0].receive_time.scale() == Time::TAI);
  REQUIRE(epoch.channels[0].transmit_time.scale() == Time::TAI);
  REQUIRE(epoch.channels[0].transmit_time < epoch.channels[0].receive_time);
  REQUIRE(epoch.channels[0].HasEphemeris());
  REQUIRE(std::isfinite(epoch.channels[0].relativistic_correction_s.val()));
  Vec3 r_sun(0.0, AU, 0.0);
  Vec3 r_rx = orbit.r();
  Vec3 r_tx = epoch.channels[0].tx_state.head(3);
  Real rx_sun_norm = (r_rx - r_sun).norm();
  Real tx_sun_norm = (r_tx - r_sun).norm();
  Real rho = (r_rx - r_tx).norm();
  Real expected_shapiro
      = 2.0 * Real(GM_SUN) / (Real(C) * Real(C))
        * log((rx_sun_norm + tx_sun_norm + rho) / (rx_sun_norm + tx_sun_norm - rho));
  REQUIRE_THAT(epoch.channels[0].shapiro_delay_m.val(), WithinAbs(expected_shapiro.val(), 1.0e-12));
  REQUIRE(epoch.values.size() == 3);
  REQUIRE(H.rows() == 3);
  REQUIRE(H.cols() == user_state.size());
  REQUIRE(epoch.covariance.rows() == 3);

  std::vector<Real> times = {0.0};
  std::vector<State> states = {user_state};
  std::vector<GNSSMeasurementsEpoch> precomputed = manager.Precompute(times, states, true);
  REQUIRE(precomputed.size() == 1);
  REQUIRE(precomputed[0].channels.size() == 1);

  std::vector<GNSSMeasurementsEpoch> precomputed_static
      = GnssMeasurement::Precompute(manager, times, states, true);
  REQUIRE(precomputed_static.size() == 1);
  REQUIRE(precomputed_static[0].channels.size() == 1);
}

TEST_CASE("measurements.gnss_measurements.moon_ci") {
  VecXd t_ephem_tai(2);
  t_ephem_tai << -100.0, 3600.0;

  MatXd rv_tx(2, 6);
  rv_tx << 10000e3, 0.0, 0.0, 0.0, 900.0, 0.0, 10000e3, 3600.0 * 900.0, 0.0, 0.0, 900.0, 0.0;

  auto constellation = MakePtr<GnssConstellation>(GnssConst::GPS);
  constellation->SetSatelliteStates({1}, t_ephem_tai, {rv_tx});

  Cart6 orbit(Vec6::Zero(), Frame::MOON_CI);
  orbit.r() = Vec3(3000e3, 0.0, 0.0);
  orbit.v() = Vec3::Zero();
  ClockState2 clock;
  JointOrbitClockState user_state(orbit, clock);

  GNSSMeasurements manager(constellation);
  GnssMeasurementOptions options;
  options.observables = {GnssObservable::PSEUDORANGE, GnssObservable::DOPPLER};
  options.receive_time_scale = Time::TDB;
  options.ephemeris_time_scale = Time::TAI;
  options.frame = Frame::MOON_CI;
  options.apply_transmitter_relativity = false;
  options.apply_shapiro_delay = false;
  options.apply_visibility = false;
  options.apply_cn0_threshold = false;
  manager.SetOptions(options);

  std::vector<Real> receive_times = {ConvertTime(Real(0.0), Time::TAI, Time::TDB)};
  std::vector<State> states = {user_state};
  std::vector<GNSSMeasurementsEpoch> precomputed = manager.Precompute(receive_times, states, true);

  REQUIRE(precomputed.size() == 1);
  REQUIRE(precomputed[0].receive_time.scale() == Time::TDB);
  REQUIRE(precomputed[0].channels.size() == 1);
  REQUIRE(precomputed[0].channels[0].frame == Frame::MOON_CI);
  REQUIRE(precomputed[0].channels[0].ephemeris_time_scale == Time::TAI);
  REQUIRE(precomputed[0].values.size() == 2);
}

TEST_CASE("measurements.gnss_measurements.batch_plasma") {
  VecXd t_ephem(2);
  t_ephem << -100.0, 3600.0;

  auto constellation = MakePtr<GnssConstellation>(GnssConst::GPS);
  constellation->SetSatelliteStates({1}, t_ephem, {CircularOrbitEci(26560e3, 0.0, t_ephem)});

  Cart6 orbit(Vec6::Zero(), Frame::ECI);
  orbit.r() = Vec3(R_EARTH + 500e3, 0.0, 0.0);
  orbit.v() = Vec3::Zero();
  ClockState2 clock;
  JointOrbitClockState user_state(orbit, clock);

  GNSSMeasurements manager(constellation);
  GnssMeasurementOptions options;
  options.apply_cn0_threshold = false;
  options.apply_ionosphere_plasma_delay = true;
  manager.SetOptions(options);
  GnssIonospherePlasmaRayTraceOptions raytrace_options;
  raytrace_options.config.straight_ray = true;
  raytrace_options.config.correction = false;
  manager.SetIonospherePlasmaRayTraceOptions(raytrace_options);
  manager.SetBatchCustomIonospherePlasmaDelayModel(
      [](const std::vector<Real>& times, const std::vector<State>& states,
         const std::vector<std::vector<GnssChannel>>& channels) {
        (void)states;
        std::vector<std::vector<Real>> delays(times.size());
        for (size_t i = 0; i < times.size(); i++) {
          delays[i].assign(channels[i].size(), 2.0 + static_cast<double>(i));
        }
        return delays;
      });

  std::vector<Real> times = {0.0, 10.0};
  std::vector<State> states = {user_state, user_state};
  std::vector<GNSSMeasurementsEpoch> precomputed = manager.Precompute(times, states, false);

  REQUIRE(precomputed.size() == 2);
  REQUIRE(precomputed[0].channels.size() == 1);
  REQUIRE(precomputed[1].channels.size() == 1);
  REQUIRE_THAT(precomputed[0].channels[0].ionosphere_plasma_delay_m.val(), WithinAbs(2.0, 1.0e-12));
  REQUIRE_THAT(precomputed[1].channels[0].ionosphere_plasma_delay_m.val(), WithinAbs(3.0, 1.0e-12));
}

TEST_CASE("measurements.gnss_measurements.visibility") {
  const Real R_earth = Real(R_EARTH);
  Vec3 r_sat(26560e3, 0.0, 0.0);

  Vec3 r_rx_visible(R_EARTH + 500e3, 0.0, 0.0);
  REQUIRE(ComputeVisibility(r_rx_visible, r_sat, R_earth));

  Vec3 r_rx_occluded(-(R_EARTH + 500e3), 0.0, 0.0);
  REQUIRE_FALSE(ComputeVisibility(r_rx_occluded, r_sat, R_earth));

  Vec3 r_surface(R_EARTH, 0.0, 0.0);
  Vec3 r_zenith_sat(26560e3, 0.0, 0.0);
  REQUIRE(ComputeVisibility(r_surface, r_zenith_sat, R_earth));

  Vec3 r_far_side_sat(-26560e3, 0.0, 0.0);
  REQUIRE_FALSE(ComputeVisibility(r_surface, r_far_side_sat, R_earth));

  Vec3 r_lunar_receiver(390000e3, 0.0, 0.0);
  Vec3 r_gps_near_side(26560e3, 0.0, 0.0);
  Vec3 r_gps_far_side(-26560e3, 0.0, 0.0);
  REQUIRE(ComputeVisibility(r_lunar_receiver, r_gps_near_side, R_earth));
  REQUIRE_FALSE(ComputeVisibility(r_lunar_receiver, r_gps_far_side, R_earth));
}
