#include <lupnt/conversions/state_conversions.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

namespace {
  // Round-trip a classical element set through Cartesian and back to Cartesian,
  // which avoids angle-wrapping ambiguity for near-singular orbits.
  void RequireCartRoundTrip(const ClassicalOE& coe, Real GM, double tol) {
    Vec6 rv = ClassicalToCart(coe, GM);
    State coe_rt = CartToClassical(Cart6(rv, coe.GetFrame()), GM);
    Vec6 rv_rt = ClassicalToCart(coe_rt, GM);
    for (int i = 0; i < 6; ++i) REQUIRE_THAT(rv_rt(i).val(), WithinAbs(rv(i).val(), tol));
  }
}  // namespace

// ---------------------------------------------------------------------------
// Near-singular orbits: near-circular and near-equatorial. The angles become
// ill-defined, but the Cartesian state must still round-trip.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.state_conversions_more.near_singular_orbits") {
  SECTION("near-circular orbit (e ~ 0)") {
    ClassicalOE coe(Vec6(7500.0e3, 1.0e-8, 45.0 * RAD, 30.0 * RAD, 60.0 * RAD, 25.0 * RAD),
                    Frame::GCRF);
    RequireCartRoundTrip(coe, GM_EARTH, 1.0e-3);
  }
  SECTION("near-equatorial orbit (i ~ 0)") {
    ClassicalOE coe(Vec6(8200.0e3, 0.05, 1.0e-8, 30.0 * RAD, 60.0 * RAD, 25.0 * RAD), Frame::GCRF);
    RequireCartRoundTrip(coe, GM_EARTH, 1.0e-3);
  }
  SECTION("near-circular AND near-equatorial") {
    ClassicalOE coe(Vec6(6900.0e3, 1.0e-9, 1.0e-9, 0.0, 0.0, 40.0 * RAD), Frame::GCRF);
    RequireCartRoundTrip(coe, GM_EARTH, 1.0e-3);
  }
}

// ---------------------------------------------------------------------------
// CartToClassical from two position vectors and a time difference (Gauss'
// method via RatioOfSectorToTriangleArea). We synthesize two positions from a
// known Keplerian orbit and confirm the recovered elements.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.state_conversions_more.cart_to_classical_two_positions") {
  Real GM = GM_EARTH;
  Real a = 9000.0e3, e = 0.1, inc = 35.0 * RAD, Om = 40.0 * RAD, w = 50.0 * RAD;
  Real M1 = 10.0 * RAD, M2 = 25.0 * RAD;

  ClassicalOE coe1(Vec6(a, e, inc, Om, w, M1), Frame::GCRF);
  ClassicalOE coe2(Vec6(a, e, inc, Om, w, M2), Frame::GCRF);

  Vec6 rv1 = ClassicalToCart(coe1, GM);
  Vec6 rv2 = ClassicalToCart(coe2, GM);

  Real n = sqrt(GM / pow(a, 3));  // mean motion
  Real dt = (M2 - M1) / n;

  Cart3 r1(Vec3(rv1.head(3)), Frame::GCRF);
  Cart3 r2(Vec3(rv2.head(3)), Frame::GCRF);

  State rec = CartToClassical(dt, r1, r2, GM);
  REQUIRE(rec.GetType() == ClassicalOE::TYPE);
  REQUIRE_THAT(rec(0).val(), WithinRel(a.val(), 1.0e-6));
  REQUIRE_THAT(rec(1).val(), WithinAbs(e.val(), 1.0e-6));
  REQUIRE_THAT(rec(2).val(), WithinAbs(inc.val(), 1.0e-6));
  REQUIRE_THAT(rec(3).val(), WithinAbs(Om.val(), 1.0e-6));
}

TEST_CASE("conversions.state_conversions_more.cart_to_classical_two_positions_equatorial") {
  // i == 0 branch: u = atan2(r1(1), r1(0)).
  Real GM = GM_EARTH;
  Real a = 8000.0e3, e = 0.05, inc = 0.0, Om = 0.0, w = 20.0 * RAD;
  Real M1 = 5.0 * RAD, M2 = 18.0 * RAD;

  ClassicalOE coe1(Vec6(a, e, inc, Om, w, M1), Frame::GCRF);
  ClassicalOE coe2(Vec6(a, e, inc, Om, w, M2), Frame::GCRF);
  Vec6 rv1 = ClassicalToCart(coe1, GM);
  Vec6 rv2 = ClassicalToCart(coe2, GM);
  Real n = sqrt(GM / pow(a, 3));
  Real dt = (M2 - M1) / n;

  State rec = CartToClassical(dt, Cart3(Vec3(rv1.head(3)), Frame::GCRF),
                              Cart3(Vec3(rv2.head(3)), Frame::GCRF), GM);
  REQUIRE_THAT(rec(0).val(), WithinRel(a.val(), 1.0e-6));
  REQUIRE_THAT(rec(1).val(), WithinAbs(e.val(), 1.0e-6));
  REQUIRE_THAT(rec(2).val(), WithinAbs(0.0, 1.0e-9));
}

// ---------------------------------------------------------------------------
// Inertial <-> synodic (chief-centered RTN) round trip.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.state_conversions_more.inertial_synodic_round_trip") {
  Cart6 rv_c(Vec3(7000.0e3, 1000.0e3, -500.0e3), Vec3(-1.2, 7.3, 0.5), Frame::GCRF);
  Cart6 rv_d(Vec3(7005.0e3, 995.0e3, -498.0e3), Vec3(-1.1, 7.31, 0.52), Frame::GCRF);

  State syn = InertialToSynodic(rv_c, rv_d);
  State back = SynodicToInertial(rv_c, syn);
  for (int i = 0; i < 6; ++i) REQUIRE_THAT(back(i).val(), WithinAbs(rv_d(i).val(), 1.0e-6));

  // The chief maps to the synodic origin (zero relative state).
  State syn_c = InertialToSynodic(rv_c, rv_c);
  for (int i = 0; i < 6; ++i) REQUIRE_THAT(syn_c(i).val(), WithinAbs(0.0, 1.0e-6));
}

// ---------------------------------------------------------------------------
// RelQuasiNonsingToClassical: absolute deputy elements from chief elements plus
// relative quasi-nonsingular differences.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.state_conversions_more.rel_quasi_nonsing_to_classical") {
  ClassicalOE coe_c(Vec6(8000.0e3, 0.02, 50.0 * RAD, 30.0 * RAD, 45.0 * RAD, 10.0 * RAD),
                    Frame::GCRF);
  // Small relative quasi-nonsingular element set [ada, adl, adex, adey, adix, adiy].
  Cart6 rqnsoe(Vec6(1000.0, 500.0, 200.0, -150.0, 0.0, 300.0), Frame::GCRF);

  State deputy = RelQuasiNonsingToClassical(coe_c, rqnsoe);
  REQUIRE(deputy.GetType() == ClassicalOE::TYPE);
  REQUIRE(deputy.GetFrame() == Frame::GCRF);
  // ad = ac + ada by construction.
  REQUIRE_THAT(deputy(0).val(), WithinAbs((coe_c(0) + rqnsoe(0)).val(), 1.0e-6));
  for (int i = 0; i < 6; ++i) REQUIRE(std::isfinite(deputy(i).val()));
}

// ---------------------------------------------------------------------------
// Vectorized (row-wise) overloads must agree with the scalar conversions.
// ---------------------------------------------------------------------------
TEST_CASE("conversions.state_conversions_more.vectorized_overloads") {
  MatX6 coe(2, 6);
  coe.row(0) << 7000.0e3, 0.05, 40.0 * RAD, 30.0 * RAD, 20.0 * RAD, 10.0 * RAD;
  coe.row(1) << 9000.0e3, 0.12, 55.0 * RAD, 10.0 * RAD, 70.0 * RAD, 80.0 * RAD;

  MatX6 rv = ClassicalToCart(coe, GM_EARTH);
  REQUIRE(rv.rows() == 2);
  for (int r = 0; r < 2; ++r) {
    Vec6 ref = ClassicalToCart(ClassicalOE(Vec6(coe.row(r).transpose()), Frame::GCRF), GM_EARTH);
    for (int c = 0; c < 6; ++c) REQUIRE_THAT(rv(r, c).val(), WithinAbs(ref(c).val(), 1.0e-3));
  }

  // Round-trip the batch back to classical elements.
  MatX6 coe_rt = CartToClassical(rv, GM_EARTH);
  REQUIRE(coe_rt.rows() == 2);
  for (int r = 0; r < 2; ++r) {
    REQUIRE_THAT(coe_rt(r, 0).val(), WithinRel(coe(r, 0).val(), 1.0e-6));  // a
    REQUIRE_THAT(coe_rt(r, 1).val(), WithinAbs(coe(r, 1).val(), 1.0e-6));  // e
  }
}
