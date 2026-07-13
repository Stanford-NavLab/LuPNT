#include <lupnt/conversions/frame_converter_spice.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <vector>

#include "../utils.cc"
#include "lupnt/interfaces/spice.h"

using namespace lupnt;
using Catch::Approx;

// Deepens frame_converter_spice.cc coverage beyond
// test_frame_converter_spice_extra.cc: additional SPICE-backed frame pairs
// (ITRF <-> MOON_CI / MOON_PA, GCRF <-> MOON_PA), the identity short-circuit,
// norm-preserving pure rotations (GCRF<->ITRF, MOON_CI<->MOON_ME), and the
// single-epoch MatX3 + multi-epoch (VecX-times) MatX6 / MatX3 / Vec3 overloads
// with their size-mismatch throws. Only Earth/Moon-relative frames are used
// (ICRF/EMR/MOON_OP need SSB-/EMB-/Sun-centered segments not in the bundled
// position cache). MOON_ME round trips are avoided: the source applies the
// same fixed B_M matrix both ways, so MOON_PA<->MOON_ME is not a true inverse.

namespace {
  Real Epoch() {
    spice::LoadSpiceKernel();
    return spice::StringToTdb("2023-04-15 00:00:00 TDB");
  }
}  // namespace

TEST_CASE("conversions.frame_converter_spice_more.round_trip_new_pairs") {
  Real et = Epoch();
  Vec6 rv;
  rv << 6800.0, -900.0, 2100.0, 0.7, 6.9, -1.1;  // km, km/s

  struct Pair {
    Frame a, b;
  };
  std::vector<Pair> pairs = {
      {Frame::ITRF, Frame::MOON_CI},
      {Frame::ITRF, Frame::MOON_PA},
      {Frame::GCRF, Frame::MOON_PA},
  };

  for (const auto& p : pairs) {
    Vec6 fwd = spice::ConvertFrameSpice(et, rv, p.a, p.b);
    Vec6 back = spice::ConvertFrameSpice(et, fwd, p.b, p.a);
    INFO("round-trip failed for a new frame pair");
    REQUIRE(back.isApprox(rv, 1e-6));
  }
}

TEST_CASE("conversions.frame_converter_spice_more.identity_conversion") {
  Real et = Epoch();
  Vec6 rv;
  rv << 5000.0, 1000.0, -2000.0, 1.0, 2.0, 3.0;

  // frame_in == frame_out short-circuits and returns the input unchanged.
  for (Frame f : {Frame::GCRF, Frame::ITRF, Frame::MOON_CI, Frame::MOON_PA}) {
    Vec6 out = spice::ConvertFrameSpice(et, rv, f, f);
    REQUIRE(out.isApprox(rv, 1e-12));
  }
}

TEST_CASE("conversions.frame_converter_spice_more.norm_preserving_rotations") {
  Real et = Epoch();
  Vec3 r(4200.0, -1700.0, 2600.0);

  // GCRF <-> ITRF is a pure rotation about the Earth center: |r| is invariant.
  Vec3 r_itrf = spice::ConvertFrameSpice(et, r, Frame::GCRF, Frame::ITRF);
  REQUIRE(r_itrf.norm().val() == Approx(r.norm().val()).epsilon(1e-9));

  // MOON_CI <-> MOON_ME is a composition of pure rotations about the Moon
  // center (MOON_CI->MOON_PA sxform, then MOON_PA->MOON_ME fixed B_M).
  Vec3 r_me = spice::ConvertFrameSpice(et, r, Frame::MOON_CI, Frame::MOON_ME);
  REQUIRE(r_me.norm().val() == Approx(r.norm().val()).epsilon(1e-9));
}

TEST_CASE("conversions.frame_converter_spice_more.matx3_single_epoch_overload") {
  Real et = Epoch();

  // Single-epoch MatX3 (many positions) overload must match per-row Vec3 calls.
  MatX3 r_in(3, 3);
  r_in << 6000.0, 0.0, 0.0, 0.0, 6000.0, 0.0, 1000.0, -2000.0, 3000.0;
  MatX3 out = spice::ConvertFrameSpice(et, r_in, Frame::GCRF, Frame::MOON_CI);
  REQUIRE(out.rows() == 3);
  for (int i = 0; i < 3; i++) {
    Vec3 expect
        = spice::ConvertFrameSpice(et, r_in.row(i).transpose().eval(), Frame::GCRF, Frame::MOON_CI);
    REQUIRE(out.row(i).transpose().isApprox(expect, 1e-9));
  }
}

TEST_CASE("conversions.frame_converter_spice_more.vecx_times_matrix_overloads") {
  Real et = Epoch();
  VecX times(3);
  times << et, et + 300.0, et + 900.0;

  // Per-epoch MatX6 states overload matches the per-row scalar conversion.
  MatX6 rv_in(3, 6);
  for (int i = 0; i < 3; i++) {
    Vec6 rv;
    rv << 7000.0 + 10.0 * i, 100.0, -300.0, 0.2, 7.4, -0.6;
    rv_in.row(i) = rv.transpose();
  }
  MatX6 out6 = spice::ConvertFrameSpice(times, rv_in, Frame::GCRF, Frame::MOON_CI);
  REQUIRE(out6.rows() == 3);
  for (int i = 0; i < 3; i++) {
    Vec6 expect = spice::ConvertFrameSpice(times(i), rv_in.row(i).transpose().eval(), Frame::GCRF,
                                           Frame::MOON_CI);
    REQUIRE(out6.row(i).transpose().isApprox(expect, 1e-9));
  }

  // Per-epoch MatX3 positions overload matches the per-row scalar Vec3 call.
  MatX3 r_in(3, 3);
  r_in << 6100.0, 0.0, 0.0, 0.0, 6200.0, 0.0, 500.0, -700.0, 900.0;
  MatX3 out3 = spice::ConvertFrameSpice(times, r_in, Frame::GCRF, Frame::MOON_CI);
  REQUIRE(out3.rows() == 3);
  for (int i = 0; i < 3; i++) {
    Vec3 expect = spice::ConvertFrameSpice(times(i), r_in.row(i).transpose().eval(), Frame::GCRF,
                                           Frame::MOON_CI);
    REQUIRE(out3.row(i).transpose().isApprox(expect, 1e-9));
  }

  // Size-mismatch guards throw.
  VecX times_bad(2);
  times_bad << et, et + 300.0;
  REQUIRE_THROWS(spice::ConvertFrameSpice(times_bad, rv_in, Frame::GCRF, Frame::MOON_CI));
  REQUIRE_THROWS(spice::ConvertFrameSpice(times_bad, r_in, Frame::GCRF, Frame::MOON_CI));
}
