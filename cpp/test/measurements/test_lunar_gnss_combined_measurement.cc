#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  GnssChannel MakeChannel(int prn, double sig_pr, double sig_dop, double sig_cp) {
    GnssChannel channel;
    channel.gnss_const = GnssConst::GPS;
    channel.prn = prn;
    channel.frequency = GnssFreq::L1;
    channel.frame = Frame::ECI;
    // Distinct transmitter positions so the two channels give different observables.
    channel.tx_state << 100.0 * prn, -50.0 * prn, 20.0 * prn, 0.0, 100.0, 0.0;
    channel.tx_clock_bias_s = 1.0e-6;
    channel.tx_clock_drift = 2.0e-9;
    channel.integer_ambiguity_cycles = 5.0;
    channel.phase_bias_cycles = 0.1;
    channel.sigma_pseudorange_m = sig_pr;
    channel.sigma_doppler_hz = sig_dop;
    channel.sigma_carrier_phase_cycles = sig_cp;
    return channel;
  }

  JointOrbitClockState MakeUserState() {
    Cart6 orbit(Vec6::Zero(), Frame::ECI);
    orbit.r() = Vec3(3.0e3, 4.0e3, 0.0);
    orbit.v() = Vec3(0.6, 0.8, 0.0);
    ClockState2 clock;
    clock.b() = 1.0e-6;
    clock.d() = 2.0e-9;
    return JointOrbitClockState(orbit, clock);
  }
}  // namespace

TEST_CASE("measurements.lunar_gnss_combined.make_tdcp_pairs") {
  std::vector<GnssChannel> current
      = {MakeChannel(1, 2.0, 0.5, 0.01), MakeChannel(2, 2.0, 0.5, 0.01)};
  std::vector<GnssChannel> previous
      = {MakeChannel(2, 2.0, 0.5, 0.01), MakeChannel(3, 2.0, 0.5, 0.01)};

  auto pairs = LunarGnssCombinedMeasurement::MakeTdcpPairs(current, previous);
  REQUIRE(pairs.size() == 1);  // only PRN 2 appears in both epochs
  REQUIRE(pairs[0].current.prn == 2);
  REQUIRE(pairs[0].previous.prn == 2);

  // A frequency mismatch breaks the match.
  std::vector<GnssChannel> prev_l5 = {MakeChannel(1, 2.0, 0.5, 0.01)};
  prev_l5[0].frequency = GnssFreq::L5;
  auto none
      = LunarGnssCombinedMeasurement::MakeTdcpPairs({MakeChannel(1, 2.0, 0.5, 0.01)}, prev_l5);
  REQUIRE(none.empty());
}

TEST_CASE("measurements.lunar_gnss_combined.measurement_covariance") {
  std::vector<GnssChannel> channels
      = {MakeChannel(1, 2.0, 0.5, 0.01), MakeChannel(2, 3.0, 0.25, 0.02)};
  GnssMeasurementOptions options;  // default observables: PR, Doppler, carrier-phase

  MatXd R = LunarGnssCombinedMeasurement::MeasurementCovariance(channels, options);
  const int n_obs = static_cast<int>(options.observables.size());
  REQUIRE(R.rows() == n_obs * 2);
  REQUIRE(R.cols() == n_obs * 2);

  // Diagonal must hold each channel's per-observable sigma^2, off-diagonal zero.
  REQUIRE_THAT(R(0, 0), WithinAbs(4.0, 1.0e-12));     // ch0 pseudorange 2.0^2
  REQUIRE_THAT(R(1, 1), WithinAbs(0.25, 1.0e-12));    // ch0 doppler 0.5^2
  REQUIRE_THAT(R(2, 2), WithinAbs(1.0e-4, 1.0e-14));  // ch0 carrier 0.01^2
  REQUIRE_THAT(R(3, 3), WithinAbs(9.0, 1.0e-12));     // ch1 pseudorange 3.0^2
  REQUIRE_THAT(R(4, 4), WithinAbs(0.0625, 1.0e-12));  // ch1 doppler 0.25^2
  REQUIRE_THAT(R(5, 5), WithinAbs(4.0e-4, 1.0e-14));  // ch1 carrier 0.02^2
  REQUIRE_THAT((R - MatXd(R.diagonal().asDiagonal())).norm(), WithinAbs(0.0, 1.0e-14));
}

TEST_CASE("measurements.lunar_gnss_combined.tdcp_covariance") {
  std::vector<GnssChannel> current = {MakeChannel(1, 2.0, 0.5, 0.01)};
  std::vector<GnssChannel> previous = {MakeChannel(1, 2.0, 0.5, 0.03)};
  auto pairs = LunarGnssCombinedMeasurement::MakeTdcpPairs(current, previous);
  REQUIRE(pairs.size() == 1);

  const double tdcp_sigma_m = 0.05;
  const double inflation_m = 0.2;
  MatXd R = LunarGnssCombinedMeasurement::TdcpCovariance(pairs, tdcp_sigma_m, inflation_m);
  REQUIRE(R.rows() == 1);

  double lambda = pairs[0].current.Wavelength().val();
  double curr_sig_m = 0.01 * lambda;
  double prev_sig_m = 0.03 * lambda;
  double expected = curr_sig_m * curr_sig_m + prev_sig_m * prev_sig_m + tdcp_sigma_m * tdcp_sigma_m
                    + inflation_m * inflation_m;
  REQUIRE_THAT(R(0, 0), WithinAbs(expected, 1.0e-12));

  // Zero inflation on the truth side is strictly smaller.
  MatXd R_truth = LunarGnssCombinedMeasurement::TdcpCovariance(pairs, tdcp_sigma_m, 0.0);
  REQUIRE(R_truth(0, 0) < R(0, 0));
}

TEST_CASE("measurements.lunar_gnss_combined.measurement_vector_matches_per_channel") {
  std::vector<GnssChannel> channels
      = {MakeChannel(1, 2.0, 0.5, 0.01), MakeChannel(2, 3.0, 0.25, 0.02)};
  GnssMeasurementOptions options;
  JointOrbitClockState x = MakeUserState();

  MatXd H;
  VecXd y = LunarGnssCombinedMeasurement::ComputeMeasurementVector(x, channels, options, &H);

  const int n_obs = static_cast<int>(options.observables.size());
  REQUIRE(y.size() == n_obs * 2);
  REQUIRE(H.rows() == n_obs * 2);
  REQUIRE(H.cols() == x.size());

  // Compare against the per-channel single measurement stacked directly.
  for (int i = 0; i < 2; ++i) {
    GnssMeasurement m(channels[i]);
    VecXd y_i = m.ComputeVector(x, nullptr, options);
    for (int j = 0; j < n_obs; ++j) REQUIRE_THAT(y(i * n_obs + j), WithinAbs(y_i(j), 1.0e-9));
  }
}

TEST_CASE("measurements.lunar_gnss_combined.compute_without_tdcp") {
  std::vector<GnssChannel> channels
      = {MakeChannel(1, 2.0, 0.5, 0.01), MakeChannel(2, 3.0, 0.25, 0.02)};
  LunarGnssCombinedMeasurement::Config cfg;
  cfg.channels = channels;
  cfg.use_tdcp = false;
  LunarGnssCombinedMeasurement meas(cfg);

  JointOrbitClockState x = MakeUserState();
  MatXd H;
  MeasData md = meas.Compute(x, &H);

  const int n_obs = static_cast<int>(cfg.current_options.observables.size());
  REQUIRE(md.value.size() == n_obs * 2);
  REQUIRE(md.covariance.rows() == n_obs * 2);
  REQUIRE(H.rows() == n_obs * 2);
  REQUIRE(H.cols() == x.size());

  // Value equals the plain stacked measurement vector.
  VecXd y_ref
      = LunarGnssCombinedMeasurement::ComputeMeasurementVector(x, channels, cfg.current_options);
  REQUIRE_THAT((md.value - y_ref).norm(), WithinAbs(0.0, 1.0e-9));

  // Covariance equals MeasurementCovariance.
  MatXd R_ref = LunarGnssCombinedMeasurement::MeasurementCovariance(channels, cfg.current_options);
  REQUIRE_THAT((md.covariance - R_ref).norm(), WithinAbs(0.0, 1.0e-12));
}

TEST_CASE("measurements.lunar_gnss_combined.compute_with_tdcp") {
  std::vector<GnssChannel> channels
      = {MakeChannel(1, 2.0, 0.5, 0.01), MakeChannel(2, 3.0, 0.25, 0.02)};
  LunarGnssCombinedMeasurement::Config cfg;
  cfg.channels = channels;
  cfg.use_tdcp = true;
  cfg.carrier_options.observables = {GnssObservable::CARRIER_PHASE};
  cfg.tdcp_pairs = LunarGnssCombinedMeasurement::MakeTdcpPairs(channels, channels);
  cfg.tdcp_sigma_m = 0.05;
  cfg.filter_tdcp_noise_inflation_m = 0.2;
  REQUIRE(cfg.tdcp_pairs.size() == 2);
  LunarGnssCombinedMeasurement meas(cfg);

  // Cloned state [x_current; x_previous].
  JointOrbitClockState x_curr = MakeUserState();
  JointOrbitClockState x_prev = MakeUserState();
  x_prev.r() = x_prev.r() + Vec3(1.0, 0.0, 0.0);  // small displacement -> nonzero TDCP
  const int base_n = static_cast<int>(x_curr.size());

  State x(2 * base_n);
  x.head(base_n) = x_curr;
  x.tail(base_n) = x_prev;

  MatXd H;
  MeasData md = meas.Compute(x, &H);

  const int n_obs = static_cast<int>(cfg.current_options.observables.size());
  const int n_pairs = static_cast<int>(cfg.tdcp_pairs.size());
  const int n_current = n_obs * 2;
  REQUIRE(md.value.size() == n_current + n_pairs);
  REQUIRE(md.covariance.rows() == n_current + n_pairs);
  REQUIRE(H.rows() == n_current + n_pairs);
  REQUIRE(H.cols() == 2 * base_n);

  // The current block only depends on the current state (right half of H is zero there).
  REQUIRE_THAT(H.block(0, base_n, n_current, base_n).norm(), WithinAbs(0.0, 1.0e-12));

  // Covariance is block-diagonal: current R then TDCP R.
  MatXd R_current
      = LunarGnssCombinedMeasurement::MeasurementCovariance(channels, cfg.current_options);
  MatXd R_tdcp = LunarGnssCombinedMeasurement::TdcpCovariance(cfg.tdcp_pairs, cfg.tdcp_sigma_m,
                                                              cfg.filter_tdcp_noise_inflation_m);
  REQUIRE_THAT((md.covariance.topLeftCorner(n_current, n_current) - R_current).norm(),
               WithinAbs(0.0, 1.0e-12));
  REQUIRE_THAT((md.covariance.bottomRightCorner(n_pairs, n_pairs) - R_tdcp).norm(),
               WithinAbs(0.0, 1.0e-12));
  // Off-diagonal cross-block is zero.
  REQUIRE_THAT(md.covariance.topRightCorner(n_current, n_pairs).norm(), WithinAbs(0.0, 1.0e-14));

  // The TDCP rows are the current-minus-previous carrier range and should be near the +x
  // displacement magnitude scale (nonzero given x_prev != x_curr).
  VecXd y_tdcp = md.value.tail(n_pairs);
  REQUIRE(y_tdcp.norm() > 0.0);
}
