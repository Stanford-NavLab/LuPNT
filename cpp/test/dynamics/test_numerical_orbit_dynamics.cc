#include <lupnt/conversions/frame_converter.h>
#include <lupnt/data/kernels.h>
#include <lupnt/dynamics/numerical_orbit_dynamics.h>
#include <lupnt/environment/forces.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;
using namespace Catch::Matchers;

const double epsilon = 1e-6;

TEST_CASE("dynamics.numerical_orbit_dynamics") {
  SECTION("NBodyDynamics scales SPICE ephemerides into configured kilometer units") {
    Real t = 0.0;
    Real t_tdb = GetLupntEpoch() + t;
    Real radius_m = 7000.0e3;
    Real speed_m_s = sqrt(GM_EARTH / radius_m);

    Vec6 rv_earth_si = GetBodyPosVel(t_tdb, BodyId::EARTH, Frame::GCRF);
    Cart6 state_si(rv_earth_si.head(3) + Vec3(radius_m, 0.0, 0.0),
                   rv_earth_si.tail(3) + Vec3(0.0, speed_m_s, 0.0), Frame::GCRF);

    NBodyDynamics dyn_si;
    dyn_si.SetFrame(Frame::GCRF);
    dyn_si.SetUseRelativity(false);
    dyn_si.AddBody(Body::Earth());

    Vec6 rv_earth_km;
    rv_earth_km.head(3) = rv_earth_si.head(3) / 1000.0;
    rv_earth_km.tail(3) = rv_earth_si.tail(3) / 1000.0;
    Cart6 state_km(rv_earth_km.head(3) + Vec3(radius_m / 1000.0, 0.0, 0.0),
                   rv_earth_km.tail(3) + Vec3(0.0, speed_m_s / 1000.0, 0.0), Frame::GCRF);

    NBodyDynamics dyn_km;
    dyn_km.SetUnits(KM_S_KG_UNITS);
    dyn_km.SetFrame(Frame::GCRF);
    dyn_km.SetUseRelativity(false);
    dyn_km.AddBody(Body::Earth(KM_S_KG_UNITS));

    Vec6 rates_si = dyn_si.ComputeRates(t, state_si);
    Vec6 rates_km = dyn_km.ComputeRates(t, state_km);

    for (int i = 0; i < 6; ++i) {
      REQUIRE_THAT(rates_km(i).val(), WithinRel(rates_si(i).val() / 1000.0, 1.0e-9));
    }
  }

  SECTION("NBodyDynamics rejects bodies in a different unit system") {
    NBodyDynamics dyn;
    dyn.SetUnits(KM_S_KG_UNITS);
    REQUIRE_THROWS(dyn.AddBody(Body::Earth()));
  }

  SECTION("JToCartTwoBodyDynamics evaluates J2 about the body-fixed spin axis") {
    Vec3 r(7000.0e3, 1200.0e3, 900.0e3);
    Cart6 state(r, Vec3::Zero(), Frame::GCRF);

    JToCartTwoBodyDynamics dyn(GM_EARTH, J2_EARTH, R_EARTH);
    VecX rates = dyn.ComputeRates(0.0, state);
    Vec3 a_j2 = rates.tail(3) + GM_EARTH * r / pow(r.norm(), 3);

    auto [R_gcrf_to_itrf, translation]
        = GetFrameRotationTranslation(GetLupntEpoch(), Frame::GCRF, Frame::ITRF);
    (void)translation;
    Vec3 r_bf = R_gcrf_to_itrf * r;
    Real r_norm = r.norm();
    Real aux1 = -3.0 / 2.0 * GM_EARTH * J2_EARTH * pow(R_EARTH, 2.0) / pow(r_norm, 5.0);
    Real aux2 = 5.0 * pow(r_bf(2) / r_norm, 2.0);
    Vec3 expected_bf(aux1 * (1.0 - aux2) * r_bf(0), aux1 * (1.0 - aux2) * r_bf(1),
                     aux1 * (3.0 - aux2) * r_bf(2));
    Vec3 expected = R_gcrf_to_itrf.transpose() * expected_bf;

    for (int i = 0; i < 3; ++i) {
      REQUIRE_THAT(a_j2(i).val(), WithinAbs(expected(i).val(), 1.0e-15));
    }
  }

  SECTION("JToCartTwoBodyDynamics can reproduce inertial-axis J2 when configured") {
    Vec3 r(7000.0e3, 1200.0e3, 900.0e3);
    Cart6 state(r, Vec3::Zero(), Frame::GCRF);

    JToCartTwoBodyDynamics dyn(GM_EARTH, J2_EARTH, R_EARTH, Frame::GCRF, Frame::GCRF);
    VecX rates = dyn.ComputeRates(0.0, state);
    Vec3 a_j2 = rates.tail(3) + GM_EARTH * r / pow(r.norm(), 3);

    Real r_norm = r.norm();
    Real aux1 = -3.0 / 2.0 * GM_EARTH * J2_EARTH * pow(R_EARTH, 2.0) / pow(r_norm, 5.0);
    Real aux2 = 5.0 * pow(r(2) / r_norm, 2.0);
    Vec3 expected(aux1 * (1.0 - aux2) * r(0), aux1 * (1.0 - aux2) * r(1),
                  aux1 * (3.0 - aux2) * r(2));

    for (int i = 0; i < 3; ++i) {
      REQUIRE_THAT(a_j2(i).val(), WithinAbs(expected(i).val(), 1.0e-15));
    }
  }

  SECTION(
      "NBodyDynamics applies relativistic correction from the Sun and closest configured body") {
    Real t = 0.0;
    Real t_tdb = GetLupntEpoch() + t;
    Real radius = 7000.0e3;
    Real speed = sqrt(GM_EARTH / radius);
    Vec3 r_rel(radius, 0.0, 0.0);
    Vec3 v_rel(0.0, speed, 0.0);

    Vec6 rv_earth = GetBodyPosVel(t_tdb, BodyId::EARTH, Frame::GCRF);
    Cart6 state(rv_earth.head(3) + r_rel, rv_earth.tail(3) + v_rel, Frame::GCRF);

    NBodyDynamics dyn_no_rel;
    dyn_no_rel.SetFrame(Frame::GCRF);
    dyn_no_rel.AddBody(Body::Earth());
    dyn_no_rel.AddBody(Body::Moon());
    dyn_no_rel.SetUseRelativity(false);

    NBodyDynamics dyn_rel;
    dyn_rel.SetFrame(Frame::GCRF);
    dyn_rel.AddBody(Body::Earth());
    dyn_rel.AddBody(Body::Moon());
    dyn_rel.SetUseRelativity(true);

    Vec6 rates_no_rel = dyn_no_rel.ComputeRates(t, state);
    Vec6 rates_rel = dyn_rel.ComputeRates(t, state);
    Vec6 rv_sun = GetBodyPosVel(t_tdb, BodyId::SUN, Frame::GCRF);
    Vec3 expected = AccelerationRelativisticCorrection(r_rel, v_rel, GM_EARTH)
                    + AccelerationRelativisticCorrection(state.head(3) - rv_sun.head(3),
                                                         state.tail(3) - rv_sun.tail(3), GM_SUN);
    Vec3 actual = rates_rel.tail(3) - rates_no_rel.tail(3);

    for (int i = 0; i < 3; ++i) {
      REQUIRE_THAT(actual(i).val(), WithinAbs(expected(i).val(), 1.0e-15));
    }
  }
}
