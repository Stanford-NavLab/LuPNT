#include <lupnt/lupnt.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

TEST_CASE("measurements.gnss_measurement_extra.channel_wavelength") {
  // Wavelength() = c / f; must match the frequency table exactly and shorten
  // monotonically with increasing carrier frequency.
  for (GnssFreq f : {GnssFreq::L1, GnssFreq::L2, GnssFreq::L5}) {
    GnssChannel ch;
    ch.frequency = f;
    Real expected = Real(C) / GNSS_FREQ_MAP.at(f);
    REQUIRE_THAT(ch.Wavelength().val(), WithinRel(expected.val(), 1e-12));
  }

  GnssChannel l1, l2, l5;
  l1.frequency = GnssFreq::L1;  // highest frequency
  l2.frequency = GnssFreq::L2;
  l5.frequency = GnssFreq::L5;  // lowest frequency
  // Higher frequency -> shorter wavelength.
  REQUIRE(l1.Wavelength().val() < l2.Wavelength().val());
  REQUIRE(l2.Wavelength().val() < l5.Wavelength().val());
  // L1 wavelength is the well-known ~0.1903 m.
  REQUIRE_THAT(l1.Wavelength().val(), WithinAbs(0.19029, 1e-4));
}

TEST_CASE("measurements.gnss_measurement_extra.channel_ephemeris") {
  SECTION("HasEphemeris false for a bare tx_state snapshot") {
    GnssChannel ch;
    ch.tx_state << 1.0e7, 0.0, 0.0, 0.0, 3.0e3, 0.0;
    REQUIRE_FALSE(ch.HasEphemeris());
    // GetTransmitState falls back to the snapshot when there is no ephemeris.
    Vec6 s = ch.GetTransmitState(Epoch::FromSeconds(Real(123.0), ch.transmit_time.scale()));
    for (int i = 0; i < 6; ++i) REQUIRE_THAT(s(i).val(), WithinAbs(ch.tx_state(i).val(), 1e-9));
  }

  SECTION("sampled state history is linearly interpolated") {
    GnssChannel ch;
    ch.transmit_time = Epoch::FromSeconds(
        Real(0.0), Time::TAI);  // keep query epoch away from the tx_state shortcut
    VecXd times(2);
    times << 0.0, 100.0;
    MatXd states(2, 6);
    states.row(0) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    states.row(1) << 100.0, 200.0, -50.0, 1.0, 2.0, 3.0;
    ch.ephemeris_times = times;
    ch.ephemeris_tx_states = states;

    REQUIRE(ch.HasEphemeris());

    // Midpoint: exactly halfway between the two samples.
    Vec6 mid = ch.GetTransmitState(Epoch::FromSeconds(Real(50.0), ch.transmit_time.scale()));
    REQUIRE_THAT(mid(0).val(), WithinAbs(50.0, 1e-9));
    REQUIRE_THAT(mid(1).val(), WithinAbs(100.0, 1e-9));
    REQUIRE_THAT(mid(2).val(), WithinAbs(-25.0, 1e-9));
    REQUIRE_THAT(mid(4).val(), WithinAbs(1.0, 1e-9));

    // Endpoints recover the samples.
    Vec6 end = ch.GetTransmitState(Epoch::FromSeconds(Real(100.0), ch.transmit_time.scale()));
    for (int i = 0; i < 6; ++i) REQUIRE_THAT(end(i).val(), WithinAbs(states(1, i), 1e-9));
  }

  SECTION("query exactly at transmit_time returns the stored snapshot") {
    GnssChannel ch;
    ch.transmit_time = Epoch::FromSeconds(Real(42.0), Time::TAI);
    ch.tx_state << 7.0, 8.0, 9.0, 0.1, 0.2, 0.3;
    // Even with an ephemeris present, an exact transmit_time hit short-circuits.
    VecXd times(2);
    times << 0.0, 100.0;
    MatXd states = MatXd::Ones(2, 6);
    ch.ephemeris_times = times;
    ch.ephemeris_tx_states = states;

    Vec6 s = ch.GetTransmitState(Epoch::FromSeconds(Real(42.0), ch.transmit_time.scale()));
    for (int i = 0; i < 6; ++i) REQUIRE_THAT(s(i).val(), WithinAbs(ch.tx_state(i).val(), 1e-12));
  }
}

TEST_CASE("measurements.gnss_measurement_extra.effective_transmit_clock") {
  GnssChannel ch;
  ch.tx_clock_bias_s = 2.5e-6;
  ch.relativistic_correction_s = -4.0e-9;
  ch.group_delay_s = 1.0e-9;
  ch.tx_clock_drift = 3.3e-12;

  // EffectiveTransmitClockBiasSeconds = bias + relativistic - group_delay.
  double expected_bias = 2.5e-6 + (-4.0e-9) - 1.0e-9;
  REQUIRE_THAT(ch.EffectiveTransmitClockBiasSeconds().val(), WithinAbs(expected_bias, 1e-18));
  REQUIRE_THAT(ch.EffectiveTransmitClockDrift().val(), WithinAbs(3.3e-12, 1e-24));
}

TEST_CASE("measurements.gnss_measurement_extra.value_as_vector") {
  GnssMeasurementValue v;
  v.pseudorange_m = 2.0e7;
  v.doppler_hz = -1500.0;
  v.carrier_phase_cycles = 1.05e8;
  v.carrier_integer_cycles = 12.0;

  SECTION("full default ordering [pseudorange, doppler, carrier]") {
    VecXd y = v.AsVector(
        {GnssObservable::PSEUDORANGE, GnssObservable::DOPPLER, GnssObservable::CARRIER_PHASE});
    REQUIRE(y.size() == 3);
    REQUIRE_THAT(y(0), WithinAbs(2.0e7, 1e-6));
    REQUIRE_THAT(y(1), WithinAbs(-1500.0, 1e-9));
    // Carrier phase includes the integer ambiguity.
    REQUIRE_THAT(y(2), WithinAbs(1.05e8 + 12.0, 1e-3));
  }

  SECTION("reordered / subset selection is honored") {
    VecXd y = v.AsVector({GnssObservable::CARRIER_PHASE, GnssObservable::PSEUDORANGE});
    REQUIRE(y.size() == 2);
    REQUIRE_THAT(y(0), WithinAbs(1.05e8 + 12.0, 1e-3));
    REQUIRE_THAT(y(1), WithinAbs(2.0e7, 1e-6));

    VecXd only_doppler = v.AsVector({GnssObservable::DOPPLER});
    REQUIRE(only_doppler.size() == 1);
    REQUIRE_THAT(only_doppler(0), WithinAbs(-1500.0, 1e-9));

    VecXd empty = v.AsVector({});
    REQUIRE(empty.size() == 0);
  }
}
