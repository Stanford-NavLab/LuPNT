#include <lupnt/conversions/frame_converter.h>
#include <lupnt/states/state.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <sstream>
#include <vector>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  const Vec6 kRv = (Vec6() << 4.2e6, -1.1e6, 2.7e6, 320.0, -410.0, 155.0).finished();

  // Frame pairs reachable via the analytic (bundled-data, no-SPICE) frame graph.
  // GCRF<->EME is a constant bias; the Moon frames use bundled DE Chebyshev
  // libration; EMR/MOON_OP use the bundled ephemeris. All are exercised through
  // the public ConvertFrame dispatcher rather than the low-level *To* helpers.
  struct FramePair {
    Frame a, b;
    double tol;
  };
}  // namespace

TEST_CASE("conversions.frame_converter_extra.predicates") {
  SECTION("IsPlanetFixedFrame / IsPlanetCiFrame / IsPlanetFrame classify correctly") {
    REQUIRE(IsPlanetFixedFrame(Frame::MARS_FIXED));
    REQUIRE(IsPlanetFixedFrame(Frame::MERCURY_FIXED));
    REQUIRE_FALSE(IsPlanetFixedFrame(Frame::MARS_CI));
    REQUIRE_FALSE(IsPlanetFixedFrame(Frame::GCRF));
    REQUIRE_FALSE(IsPlanetFixedFrame(Frame::MOON_PA));  // Moon is not a "planet" frame here

    REQUIRE(IsPlanetCiFrame(Frame::MARS_CI));
    REQUIRE(IsPlanetCiFrame(Frame::NEPTUNE_CI));
    REQUIRE_FALSE(IsPlanetCiFrame(Frame::MARS_FIXED));
    REQUIRE_FALSE(IsPlanetCiFrame(Frame::MOON_CI));

    REQUIRE(IsPlanetFrame(Frame::MARS_CI));
    REQUIRE(IsPlanetFrame(Frame::MARS_FIXED));
    REQUIRE_FALSE(IsPlanetFrame(Frame::GCRF));
    REQUIRE_FALSE(IsPlanetFrame(Frame::MOON_CI));
  }

  SECTION("GetFrameCenter covers planet frames") {
    REQUIRE(GetFrameCenter(Frame::MARS_CI) == BodyId::MARS);
    REQUIRE(GetFrameCenter(Frame::MARS_FIXED) == BodyId::MARS);
    REQUIRE(GetFrameCenter(Frame::JUPITER_CI) == BodyId::JUPITER);
  }
}

TEST_CASE("conversions.frame_converter_extra.stream_operator") {
  auto to_str = [](Frame f) {
    std::ostringstream os;
    os << f;
    return os.str();
  };
  REQUIRE(to_str(Frame::GCRF) == "GCRF");
  REQUIRE(to_str(Frame::MOON_PA) == "MOON_PA");
  REQUIRE(to_str(Frame::ITRF) == "ITRF");
}

TEST_CASE("conversions.frame_converter_extra.same_frame_identity") {
  std::vector<Frame> frames
      = {Frame::GCRF, Frame::EME, Frame::MOON_CI, Frame::MOON_PA, Frame::MOON_ME, Frame::MARS_CI};
  for (Frame f : frames) {
    Vec6 rv_out = ConvertFrame(Real(0.0), kRv, f, f);
    for (int i = 0; i < 6; ++i) REQUIRE_THAT(rv_out(i).val(), WithinAbs(kRv(i).val(), 1e-9));
    Vec3 r_out = ConvertFrame(Real(0.0), Vec3(kRv.head(3)), f, f);
    for (int i = 0; i < 3; ++i) REQUIRE_THAT(r_out(i).val(), WithinAbs(kRv(i).val(), 1e-9));
  }
}

TEST_CASE("conversions.frame_converter_extra.roundtrip_and_overloads") {
  std::vector<FramePair> pairs = {
      {Frame::GCRF, Frame::EME, 1e-6},        {Frame::MOON_CI, Frame::MOON_PA, 1e-6},
      {Frame::MOON_PA, Frame::MOON_ME, 1e-6}, {Frame::MOON_CI, Frame::MOON_OP, 1e-6},
      {Frame::GCRF, Frame::EMR, 1e-5},        {Frame::GCRF, Frame::MOON_CI, 1e-5},
  };

  for (const auto& p : pairs) {
    INFO("pair " << p.a << " <-> " << p.b);
    Real t(3.0 * 86400.0);

    SECTION("Vec6 A->B->A round-trips") {
      Vec6 rv_b = ConvertFrame(t, kRv, p.a, p.b);
      Vec6 rv_a = ConvertFrame(t, rv_b, p.b, p.a);
      for (int i = 0; i < 6; ++i) REQUIRE_THAT(rv_a(i).val(), WithinAbs(kRv(i).val(), p.tol));
    }

    SECTION("Vec3 and Vec6 overloads agree on the position block") {
      Vec3 r_v3 = ConvertFrame(t, Vec3(kRv.head(3)), p.a, p.b);
      Vec6 rv_v6 = ConvertFrame(t, kRv, p.a, p.b);
      for (int i = 0; i < 3; ++i) REQUIRE_THAT(rv_v6(i).val(), WithinAbs(r_v3(i).val(), p.tol));
    }
  }
}

TEST_CASE("conversions.frame_converter_extra.rotation_translation") {
  // GetFrameRotationTranslation returns (R, t) with r_to = R*r_from + t. The
  // rotation block must be orthonormal, and the affine map must reproduce the
  // ConvertFrame position result.
  auto proper_rotation = [](const Mat3& R, double tol) {
    Mat3 I = R * R.transpose();
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) REQUIRE_THAT(I(i, j).val(), WithinAbs(i == j ? 1.0 : 0.0, tol));
    REQUIRE_THAT(R.determinant().val(), WithinAbs(1.0, tol));
  };

  std::vector<FramePair> pairs = {
      {Frame::GCRF, Frame::EME, 1e-6},
      {Frame::MOON_CI, Frame::MOON_PA, 1e-6},
      {Frame::GCRF, Frame::MOON_CI, 1e-5},
  };

  for (const auto& p : pairs) {
    INFO("pair " << p.a << " -> " << p.b);
    Real t(2.0 * 86400.0);

    auto [R, tr] = GetFrameRotationTranslation(t, p.a, p.b);
    proper_rotation(R, 1e-9);

    // The affine (R, tr) must reproduce the full ConvertFrame position (the
    // affine map is derived from basis-vector transforms of the same chain).
    Vec3 r_affine = R * Vec3(kRv.head(3)) + tr;
    Vec3 r_convert = ConvertFrame(t, Vec3(kRv.head(3)), p.a, p.b);
    for (int i = 0; i < 3; ++i)
      REQUIRE_THAT(r_affine(i).val(), WithinAbs(r_convert(i).val(), p.tol));

    // The `rotate_only` position overload applies exactly this affine (R, tr)
    // from GetFrameRotationTranslation, so it must match r_affine.
    Vec3 r_rot_only = ConvertFrame(t, Vec3(kRv.head(3)), p.a, p.b, /*rotate_only=*/true);
    for (int i = 0; i < 3; ++i)
      REQUIRE_THAT(r_rot_only(i).val(), WithinAbs(r_affine(i).val(), 1e-6));
  }
}

TEST_CASE("conversions.frame_converter_extra.rotation_translation_rv") {
  // GetFrameRotationTranslationRv returns the 6x6 affine (R6, t6) with
  // x_to = R6*x_from + t6, which must reproduce the Vec6 ConvertFrame result
  // (including the velocity block's rotation-rate coupling).
  std::vector<FramePair> pairs = {
      {Frame::GCRF, Frame::EME, 1e-4},
      {Frame::MOON_CI, Frame::MOON_PA, 1e-4},
      {Frame::GCRF, Frame::MOON_CI, 1e-3},
  };

  for (const auto& p : pairs) {
    INFO("pair " << p.a << " -> " << p.b);
    Real t(4.0 * 86400.0);

    auto [R6, t6] = GetFrameRotationTranslationRv(t, p.a, p.b);
    REQUIRE(R6.rows() == 6);
    REQUIRE(R6.cols() == 6);

    Vec6 x_affine = R6 * kRv + t6;
    Vec6 x_convert = ConvertFrame(t, kRv, p.a, p.b);
    for (int i = 0; i < 6; ++i)
      REQUIRE_THAT(x_affine(i).val(), WithinAbs(x_convert(i).val(), p.tol));

    // The upper-left 3x3 block of R6 is the position rotation and must match
    // GetFrameRotationTranslation's R.
    auto [R3, tr3] = GetFrameRotationTranslation(t, p.a, p.b);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) REQUIRE_THAT(R6(i, j).val(), WithinAbs(R3(i, j).val(), 1e-9));
  }
}

TEST_CASE("conversions.frame_converter_extra.matrix_and_multiepoch") {
  Real t(1.5 * 86400.0);
  Frame in = Frame::GCRF, out = Frame::MOON_CI;

  SECTION("MatX6 row-wise overload matches per-row Vec6 conversion") {
    MatX6 rv(3, 6);
    rv.row(0) = kRv.transpose();
    rv.row(1) = (2.0 * kRv).transpose();
    rv.row(2) = (-0.5 * kRv).transpose();

    MatX6 out_mat = ConvertFrame(t, rv, in, out);
    REQUIRE(out_mat.rows() == 3);
    for (int r = 0; r < 3; ++r) {
      Vec6 ref = ConvertFrame(t, Vec6(rv.row(r).transpose()), in, out);
      for (int c = 0; c < 6; ++c) REQUIRE_THAT(out_mat(r, c).val(), WithinAbs(ref(c).val(), 1e-4));
    }
  }

  SECTION("MatX3 row-wise overload matches per-row Vec3 conversion") {
    MatX3 r(2, 3);
    r.row(0) = kRv.head(3).transpose();
    r.row(1) = (3.0 * kRv.head(3)).transpose();

    MatX3 out_mat = ConvertFrame(t, r, in, out);
    REQUIRE(out_mat.rows() == 2);
    for (int rr = 0; rr < 2; ++rr) {
      Vec3 ref = ConvertFrame(t, Vec3(r.row(rr).transpose()), in, out);
      for (int c = 0; c < 3; ++c) REQUIRE_THAT(out_mat(rr, c).val(), WithinAbs(ref(c).val(), 1e-4));
    }
  }

  SECTION("multi-epoch (VecX, single Vec6) overload converts per-epoch") {
    VecX epochs(3);
    epochs << 0.0, 2.0 * 86400.0, 5.0 * 86400.0;
    MatX6 out_mat = ConvertFrame(epochs, kRv, in, out);
    REQUIRE(out_mat.rows() == 3);
    for (int k = 0; k < 3; ++k) {
      Vec6 ref = ConvertFrame(Real(epochs(k)), kRv, in, out);
      for (int c = 0; c < 6; ++c) REQUIRE_THAT(out_mat(k, c).val(), WithinAbs(ref(c).val(), 1e-4));
    }
  }

  SECTION("time-tagged rows (VecX, MatX6) overload uses per-row epoch") {
    VecX epochs(2);
    epochs << 0.0, 3.0 * 86400.0;
    MatX6 rv(2, 6);
    rv.row(0) = kRv.transpose();
    rv.row(1) = (2.0 * kRv).transpose();

    MatX6 out_mat = ConvertFrame(epochs, rv, in, out);
    REQUIRE(out_mat.rows() == 2);
    for (int r = 0; r < 2; ++r) {
      Vec6 ref = ConvertFrame(Real(epochs(r)), Vec6(rv.row(r).transpose()), in, out);
      for (int c = 0; c < 6; ++c) REQUIRE_THAT(out_mat(r, c).val(), WithinAbs(ref(c).val(), 1e-4));
    }
  }
}

TEST_CASE("conversions.frame_converter_extra.state_overload") {
  // The typed-State ConvertFrame overload round-trips through Cartesian and back
  // to the input StateType, preserving the representation.
  Real t(2.0 * 86400.0);
  Cart6 x_in(kRv, Frame::GCRF);
  State x_moon = ConvertFrame(t, State(x_in), Frame::MOON_CI);
  State x_back = ConvertFrame(t, x_moon, Frame::GCRF);
  for (int i = 0; i < 6; ++i) REQUIRE_THAT(x_back(i).val(), WithinAbs(kRv(i).val(), 1e-4));
}
