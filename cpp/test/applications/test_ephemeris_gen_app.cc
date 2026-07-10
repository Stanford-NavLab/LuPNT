#include <lupnt/agents/satellite.h>
#include <lupnt/applications/ephemeris/ephemeris_gen_app.h>
#include <lupnt/applications/ephemeris/lunanet_sat_app.h>
#include <lupnt/dynamics/surface_dynamics.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../utils.cc"

using namespace lupnt;

// A short constant-velocity (straight-line) arc: a pure polynomial ephemeris
// (no Keplerian baseline) reproduces it to machine precision, so these tests are
// fully self-contained and need no gravity/ephemeris data files.
namespace {
  std::pair<VecXd, MatXd> StraightLineArc(int n, double dt) {
    VecXd t_s(n);
    MatXd rv(n, 6);
    Vec6d x0;
    x0 << 7.0e6, 1.0e5, -2.0e5, 10.0, -5.0, 3.0;  // [m, m, m, m/s, m/s, m/s]
    for (int i = 0; i < n; ++i) {
      const double t = dt * i;
      t_s(i) = t;
      Vec6d x = x0;
      x.head(3) += x0.tail(3) * t;
      rv.row(i) = x.transpose();
    }
    return {t_s, rv};
  }
}  // namespace

TEST_CASE("applications.ephemeris_gen_app.explicit_arc") {
  auto [t_s, rv] = StraightLineArc(61, 60.0);

  EphemerisGenConfig cfg;
  cfg.ephemeris_options.use_keplerian_baseline = false;
  cfg.ephemeris_options.poly_order = 3;
  cfg.almanac_options.poly_order = 2;
  EphemerisGenApp app(cfg);

  const BroadcastMessage& msg = app.GenerateEphemerisFromArc(0.0, t_s, rv);
  REQUIRE(app.GetEphemerisMessages().size() == 1);
  REQUIRE(msg.params.size() == app.GetEphemerisModel().NumParams());
  REQUIRE_THAT(msg.t_start_s, Catch::Matchers::WithinAbs(t_s(0), 1e-9));
  REQUIRE_THAT(msg.t_end_s, Catch::Matchers::WithinAbs(t_s(t_s.size() - 1), 1e-9));

  // The fitted ephemeris reproduces the straight-line arc.
  MatXd fit = app.GetEphemerisModel().Eval(t_s, msg.params);
  REQUIRE((fit.leftCols(3) - rv.leftCols(3)).cwiseAbs().maxCoeff() < 1e-3);

  // Validity-window (latest-page) selection.
  REQUIRE(app.LatestEphemeris(1800.0) != nullptr);
  REQUIRE(app.LatestEphemeris(t_s(t_s.size() - 1) + 1.0e6) == nullptr);

  // Almanac from the same arc appends to its own message list.
  app.GenerateAlmanacFromArc(0.0, t_s, rv);
  REQUIRE(app.GetAlmanacMessages().size() == 1);
}

TEST_CASE("applications.ephemeris_gen_app.agent_driven_step") {
  // A satellite with a constant (static) state: GetStateAt returns the same
  // Cartesian state at every epoch, so a constant/linear polynomial ephemeris
  // reproduces the sampled arc exactly.
  Satellite sat;
  sat.SetName("lunanet_sat_provider");
  sat.SetTime(0.0);
  sat.SetDynamics(MakePtr<StaticDynamics>());
  Cart6 state(Vec6(7.0e6, 1.0e5, -2.0e5, 10.0, -5.0, 3.0), Frame::MOON_CI);
  sat.SetState(state);

  EphemerisGenConfig cfg;
  cfg.generate_almanac = false;  // a static state is not a valid Kepler orbit
  cfg.ephemeris_options.use_keplerian_baseline = false;
  cfg.ephemeris_options.poly_order = 1;
  cfg.ephemeris_window_s = 3600.0;
  cfg.ephemeris_fit_samples = 11;

  auto gen = MakePtr<EphemerisGenApp>(cfg);

  // Embed the generator in the LunaNet satellite routine.
  LunaNetSatApp sat_app;
  sat_app.SetName("lunanet_sat");
  sat_app.SetAgent(&sat);
  sat_app.AddSubApp(gen);

  sat_app.Setup();
  sat_app.Step(0.0);     // generation due at t=0 -> one ephemeris message
  sat_app.Step(1800.0);  // still inside the first refresh window -> no new message
  sat_app.Finish();

  REQUIRE(gen->GetEphemerisMessages().size() == 1);
  REQUIRE(gen->GetAlmanacMessages().empty());

  const BroadcastMessage& msg = gen->GetEphemerisMessages().front();
  REQUIRE(msg.frame == Frame::MOON_CI);
  REQUIRE_THAT(msg.t_start_s, Catch::Matchers::WithinAbs(0.0, 1e-9));
  REQUIRE_THAT(msg.t_end_s, Catch::Matchers::WithinAbs(3600.0, 1e-9));

  // The generated ephemeris reproduces the (constant) broadcast state.
  VecXd t_q(3);
  t_q << 0.0, 1800.0, 3600.0;
  MatXd fit = gen->GetEphemerisModel().Eval(t_q, msg.params);
  for (int i = 0; i < t_q.size(); ++i) {
    REQUIRE((fit.row(i).head(3).transpose() - state.r().cast<double>()).cwiseAbs().maxCoeff()
            < 1e-3);
  }

  // A receiver-style lookup returns the page valid over the window.
  REQUIRE(gen->LatestEphemeris(1800.0) == &msg);
  REQUIRE(gen->LatestEphemeris(7200.0) == nullptr);
}
