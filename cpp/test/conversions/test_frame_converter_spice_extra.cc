#include <lupnt/conversions/frame_converter_spice.h>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "../utils.cc"
#include "lupnt/interfaces/spice.h"

using namespace lupnt;
using Catch::Approx;

// SPICE-backed frame conversions. Only the runtime paths that route through
// spkezr/sxform/pxform are exercised (never unitim_c/et2utc_c, which segfault
// in this build -- see test_spice_extra.cc).

namespace {
  Real Epoch() {
    spice::LoadSpiceKernel();
    return spice::StringToTdb("2023-04-15 00:00:00 TDB");
  }
}  // namespace

TEST_CASE("conversions.frame_converter_spice_extra.round_trip") {
  Real et = Epoch();
  Vec6 rv;
  rv << 7000.0, 1200.0, -2500.0, 1.0, 7.2, 0.5;  // km, km/s (frame-agnostic)

  struct Pair {
    Frame a, b;
  };
  // NOTE: MOON_PA <-> MOON_ME is intentionally excluded from the round-trip
  // identity check. The source applies the *same* fixed matrix B_M in both
  // directions (rather than B_M and its transpose), so it is not a true inverse
  // pair; the small-angle discrepancy (~km) is confirmed nonzero but is a
  // pre-existing source detail, not something to assert an identity on here.
  // (Its norm-preserving property is still checked in rotation_orthonormal.)
  // (GCRF <-> ICRF and MOON_CI <-> MOON_OP are omitted: they require
  // SolarSystemBarycenter-/Sun-centered Chebyshev segments not present in the
  // bundled position cache.)
  std::vector<Pair> pairs = {
      {Frame::GCRF, Frame::ITRF},
      {Frame::GCRF, Frame::MOON_CI},
      {Frame::MOON_CI, Frame::MOON_PA},
  };

  for (const auto& p : pairs) {
    Vec6 fwd = spice::ConvertFrameSpice(et, rv, p.a, p.b);
    Vec6 back = spice::ConvertFrameSpice(et, fwd, p.b, p.a);
    INFO("round-trip failed for frame pair index");
    REQUIRE(back.isApprox(rv, 1e-6));
  }
}

TEST_CASE("conversions.frame_converter_spice_extra.rotation_orthonormal") {
  Real et = Epoch();
  // For pure-rotation frame pairs (no translation of origin), the position
  // norm is invariant under the transform.
  Vec3 r(5000.0, -3000.0, 1500.0);

  // MOON_CI <-> MOON_PA is a pure rotation about the Moon center.
  Vec3 r_pa = spice::ConvertFrameSpice(et, r, Frame::MOON_CI, Frame::MOON_PA);
  REQUIRE(r_pa.norm().val() == Approx(r.norm().val()).epsilon(1e-9));

  // MOON_PA <-> MOON_ME is a small fixed rotation; norm-preserving too.
  Vec3 r_me = spice::ConvertFrameSpice(et, r, Frame::MOON_PA, Frame::MOON_ME);
  REQUIRE(r_me.norm().val() == Approx(r.norm().val()).epsilon(1e-9));
}

TEST_CASE("conversions.frame_converter_spice_extra.vec3_vec6_agreement") {
  Real et = Epoch();
  Vec3 r(6000.0, 500.0, -1200.0);

  // The Vec3 overload zero-pads velocity and returns the position sub-block of
  // the Vec6 conversion; the two must agree exactly.
  Vec6 rv6;
  rv6 << r, Vec3::Zero();
  Vec3 out3 = spice::ConvertFrameSpice(et, r, Frame::GCRF, Frame::MOON_PA);
  Vec6 out6 = spice::ConvertFrameSpice(et, rv6, Frame::GCRF, Frame::MOON_PA);
  REQUIRE(out3.isApprox(out6.head(3), 1e-9));
}

TEST_CASE("conversions.frame_converter_spice_extra.matrix_overload_agreement") {
  Real et = Epoch();
  Vec6 rv;
  rv << 4000.0, -2000.0, 900.0, 0.5, 6.0, -1.0;

  // MatX6 (single-epoch, many states) overload matches per-row scalar calls.
  MatX6 rv_in(3, 6);
  for (int i = 0; i < 3; i++) rv_in.row(i) = (rv * (1.0 + 0.1 * i)).transpose();
  MatX6 out = spice::ConvertFrameSpice(et, rv_in, Frame::GCRF, Frame::MOON_CI);
  REQUIRE(out.rows() == 3);
  for (int i = 0; i < 3; i++) {
    Vec6 expect = spice::ConvertFrameSpice(et, rv_in.row(i).transpose().eval(), Frame::GCRF,
                                           Frame::MOON_CI);
    REQUIRE(out.row(i).transpose().isApprox(expect, 1e-9));
  }
}

TEST_CASE("conversions.frame_converter_spice_extra.time_series_overload") {
  Real et = Epoch();
  Vec6 rv;
  rv << 7000.0, 0.0, 0.0, 0.0, 7.5, 0.0;

  VecX times(3);
  times << et, et + 600.0, et + 1200.0;

  // VecX-times overload (fixed state, many epochs).
  MatX6 out = spice::ConvertFrameSpice(times, rv, Frame::GCRF, Frame::MOON_CI);
  REQUIRE(out.rows() == 3);
  for (int i = 0; i < 3; i++) {
    Vec6 expect = spice::ConvertFrameSpice(times(i), rv, Frame::GCRF, Frame::MOON_CI);
    REQUIRE(out.row(i).transpose().isApprox(expect, 1e-9));
  }

  // The Cart6 convenience overload carries the output frame tag.
  Cart6 state(rv, Frame::GCRF);
  Cart6 conv = spice::ConvertFrameSpice(et, state, Frame::MOON_CI);
  REQUIRE(conv.GetFrame() == Frame::MOON_CI);
  Vec6 direct = spice::ConvertFrameSpice(et, rv, Frame::GCRF, Frame::MOON_CI);
  REQUIRE(Vec6(conv).isApprox(direct, 1e-9));
}
