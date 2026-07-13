#include <lupnt/conversions/frame_converter.h>
#include <lupnt/states/state.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  const Vec6 kRv = (Vec6() << 3.9e6, 8.0e5, -2.1e6, 210.0, -330.0, 90.0).finished();
}

// ---------------------------------------------------------------------------
// GetFrameCenter over the Earth/Moon/SSB frames (test_frame_converter_extra.cc
// only covers the planet frames), plus the not-found error path.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.frame_converter_more.frame_center_earth_moon_ssb") {
  for (Frame f : {Frame::GCRF, Frame::ITRF, Frame::EME, Frame::ECEF, Frame::ECI})
    REQUIRE(GetFrameCenter(f) == BodyId::EARTH);

  for (Frame f : {Frame::MOON_CI, Frame::MOON_PA, Frame::MOON_ME, Frame::MOON_OP})
    REQUIRE(GetFrameCenter(f) == BodyId::MOON);

  REQUIRE(GetFrameCenter(Frame::ICRF) == BodyId::SOLAR_SYSTEM_BARYCENTER);

  // A frame with no entry in frame_centers triggers the LUPNT_CHECK failure.
  REQUIRE_THROWS(GetFrameCenter(Frame::UNDEFINED));
}

// ---------------------------------------------------------------------------
// Multi-epoch position overloads that test_frame_converter_extra.cc does not
// exercise: (VecX, Vec3) and (VecX, MatX3). Analytic GCRF<->MOON_CI path uses
// only bundled ephemeris (no SPICE).
// ---------------------------------------------------------------------------
TEST_CASE("conversions.frame_converter_more.multi_epoch_position_overloads") {
  Frame in = Frame::GCRF, out = Frame::MOON_CI;
  Vec3 r = kRv.head(3);

  SECTION("(VecX, Vec3): one converted position per epoch") {
    VecX epochs(3);
    epochs << 0.0, 2.0 * 86400.0, 5.0 * 86400.0;
    MatX3 out_mat = ConvertFrame(epochs, r, in, out);
    REQUIRE(out_mat.rows() == 3);
    REQUIRE(out_mat.cols() == 3);
    for (int k = 0; k < 3; ++k) {
      Vec3 ref = ConvertFrame(Real(epochs(k)), r, in, out);
      for (int c = 0; c < 3; ++c) REQUIRE_THAT(out_mat(k, c).val(), WithinAbs(ref(c).val(), 1e-4));
    }
  }

  SECTION("(VecX, MatX3): per-row epoch applied to per-row position") {
    VecX epochs(2);
    epochs << 0.0, 3.0 * 86400.0;
    MatX3 r_in(2, 3);
    r_in.row(0) = r.transpose();
    r_in.row(1) = (2.0 * r).transpose();

    MatX3 out_mat = ConvertFrame(epochs, r_in, in, out);
    REQUIRE(out_mat.rows() == 2);
    for (int row = 0; row < 2; ++row) {
      Vec3 ref = ConvertFrame(Real(epochs(row)), Vec3(r_in.row(row).transpose()), in, out);
      for (int c = 0; c < 3; ++c)
        REQUIRE_THAT(out_mat(row, c).val(), WithinAbs(ref(c).val(), 1e-4));
    }
  }
}

// ---------------------------------------------------------------------------
// Vec6 rotate_only path (extra only covers the Vec3 rotate_only). The Vec6
// rotate_only position block applies the same affine (R, r) as the Vec3
// overload and as the full transform; the velocity block instead uses only the
// rotation R (no frame-origin velocity or rotation-rate coupling), so it
// differs from the full transform's velocity.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.frame_converter_more.vec6_rotate_only") {
  Real t(2.5 * 86400.0);
  Frame in = Frame::GCRF, out = Frame::MOON_CI;

  Vec6 rot_only = ConvertFrame(t, kRv, in, out, /*rotate_only=*/true);
  Vec3 rot_only_v3 = ConvertFrame(t, Vec3(kRv.head(3)), in, out, /*rotate_only=*/true);
  Vec6 full = ConvertFrame(t, kRv, in, out);

  // Position block matches the Vec3 rotate_only overload exactly...
  for (int i = 0; i < 3; ++i)
    REQUIRE_THAT(rot_only(i).val(), WithinAbs(rot_only_v3(i).val(), 1e-6));
  // ...and equals the full transform's position (the affine r includes the
  // frame-origin translation).
  for (int i = 0; i < 3; ++i) REQUIRE_THAT(rot_only(i).val(), WithinAbs(full(i).val(), 1e-3));

  // The velocity block is rotation-only, so it drops the Earth->Moon relative
  // velocity (~1 km/s) present in the full transform.
  double vel_diff = (full.tail(3) - rot_only.tail(3)).norm().val();
  REQUIRE(vel_diff > 1.0);

  // Rotation preserves the input speed in the rotate-only velocity block.
  REQUIRE_THAT(rot_only.tail(3).norm().val(), WithinRel(kRv.tail(3).norm().val(), 1e-6));
}

// ---------------------------------------------------------------------------
// Same-frame identity for the affine transform getters: from == to must give
// R = I and t = 0.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.frame_converter_more.same_frame_affine_identity") {
  Real t(1.0 * 86400.0);
  for (Frame f : {Frame::GCRF, Frame::MOON_CI, Frame::EME}) {
    auto [R, tr] = GetFrameRotationTranslation(t, f, f);
    for (int i = 0; i < 3; ++i) {
      REQUIRE_THAT(tr(i).val(), WithinAbs(0.0, 1e-9));
      for (int j = 0; j < 3; ++j) REQUIRE_THAT(R(i, j).val(), WithinAbs(i == j ? 1.0 : 0.0, 1e-9));
    }

    auto [R6, t6] = GetFrameRotationTranslationRv(t, f, f);
    for (int i = 0; i < 6; ++i) {
      REQUIRE_THAT(t6(i).val(), WithinAbs(0.0, 1e-9));
      REQUIRE_THAT(R6(i, i).val(), WithinAbs(1.0, 1e-9));
    }
  }
}
