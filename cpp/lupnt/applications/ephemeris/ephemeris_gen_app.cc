#include "lupnt/applications/ephemeris/ephemeris_gen_app.h"

#include <fmt/format.h>

#include <utility>

#include "lupnt/agents/agent.h"
#include "lupnt/core/logger.h"
#include "lupnt/states/state.h"

namespace lupnt {

  namespace {
    // Small tolerance so a refresh scheduled exactly at the current step fires
    // despite floating-point round-off in the accumulated step times.
    constexpr double kTimeEps = 1.0e-6;
  }  // namespace

  EphemerisGenApp::EphemerisGenApp(EphemerisGenConfig config)
      : LunaNetSubApp("ephemeris_gen"),
        config_(std::move(config)),
        ephemeris_(config_.ephemeris_options),
        almanac_(config_.almanac_options) {}

  void EphemerisGenApp::Setup(LunaNetSatApp& app) {
    LunaNetSubApp::Setup(app);
    sat_app_ = &app;
    // Broadcast in the configured output frame: the fit models interpret the
    // (converted) arc in this frame and evaluate their output in it.
    config_.ephemeris_options.frame = config_.output_frame;
    config_.almanac_options.frame = config_.output_frame;
    ephemeris_ = CartesianEphemeris(config_.ephemeris_options);
    almanac_ = Almanac(config_.almanac_options);
    ephemeris_msgs_.clear();
    almanac_msgs_.clear();
    next_ephemeris_gen_s_ = 0.0;
    next_almanac_gen_s_ = 0.0;
  }

  void EphemerisGenApp::Step(Real t) {
    const double td = t.val();
    Agent* agent = sat_app_ != nullptr ? sat_app_->GetAgent() : nullptr;
    // Without an agent there is no predicted arc to sample; the models can still
    // be driven explicitly via Generate*FromArc.
    if (agent == nullptr) return;

    if (config_.generate_ephemeris && td + kTimeEps >= next_ephemeris_gen_s_) {
      GenerateFromAgent(td, config_.ephemeris_window_s, config_.ephemeris_fit_samples, true);
      const double refresh = config_.ephemeris_refresh_s > 0.0 ? config_.ephemeris_refresh_s
                                                               : config_.ephemeris_window_s;
      next_ephemeris_gen_s_ = td + refresh;
    }
    if (config_.generate_almanac && td + kTimeEps >= next_almanac_gen_s_) {
      GenerateFromAgent(td, config_.almanac_window_s, config_.almanac_fit_samples, false);
      const double refresh
          = config_.almanac_refresh_s > 0.0 ? config_.almanac_refresh_s : config_.almanac_window_s;
      next_almanac_gen_s_ = td + refresh;
    }
  }

  void EphemerisGenApp::Finish() {
    LunaNetSubApp::Finish();
    Logger::Debug(fmt::format("{}: generated {} ephemeris and {} almanac message(s)", GetName(),
                              ephemeris_msgs_.size(), almanac_msgs_.size()),
                  "EphemerisGenApp");
  }

  const BroadcastMessage& EphemerisGenApp::GenerateEphemerisFromArc(double t_gen_s,
                                                                    const VecXd& t_s,
                                                                    const MatXd& rv) {
    LUPNT_CHECK(t_s.size() == rv.rows() && rv.cols() == 6,
                "Arc shape mismatch: t_s and rv must have the same length and rv must be N x 6",
                "EphemerisGenApp");
    BroadcastMessage msg;
    msg.t_generated_s = t_gen_s;
    msg.t_start_s = t_s(0);
    msg.t_end_s = t_s(t_s.size() - 1);
    msg.params = ephemeris_.Fit(t_s, rv);
    ephemeris_msgs_.push_back(std::move(msg));
    return ephemeris_msgs_.back();
  }

  const BroadcastMessage& EphemerisGenApp::GenerateAlmanacFromArc(double t_gen_s, const VecXd& t_s,
                                                                  const MatXd& rv) {
    LUPNT_CHECK(t_s.size() == rv.rows() && rv.cols() == 6,
                "Arc shape mismatch: t_s and rv must have the same length and rv must be N x 6",
                "EphemerisGenApp");
    BroadcastMessage msg;
    msg.t_generated_s = t_gen_s;
    msg.t_start_s = t_s(0);
    msg.t_end_s = t_s(t_s.size() - 1);
    msg.params = almanac_.Fit(t_s, rv);
    almanac_msgs_.push_back(std::move(msg));
    return almanac_msgs_.back();
  }

  const BroadcastMessage* EphemerisGenApp::LatestEphemeris(double t_s) const {
    return Latest(ephemeris_msgs_, t_s);
  }

  const BroadcastMessage* EphemerisGenApp::LatestAlmanac(double t_s) const {
    return Latest(almanac_msgs_, t_s);
  }

  void EphemerisGenApp::GenerateFromAgent(double t_gen_s, double window_s, int n_samples,
                                          bool is_ephemeris) {
    Agent* agent = sat_app_ != nullptr ? sat_app_->GetAgent() : nullptr;
    LUPNT_CHECK(agent != nullptr,
                "EphemerisGenApp requires an owning agent to sample the predicted arc",
                "EphemerisGenApp");
    LUPNT_CHECK(n_samples >= 2, "fit_samples must be >= 2", "EphemerisGenApp");
    LUPNT_CHECK(window_s > 0.0, "validity window must be positive", "EphemerisGenApp");

    // Sample the agent's predicted state over the validity window, i.e. propagate
    // the satellite's own orbit forward the way a real vehicle would before
    // broadcasting an ephemeris/almanac for the upcoming window.
    VecXd t_s(n_samples);
    MatXd rv(n_samples, 6);
    Frame agent_frame = Frame::MOON_CI;
    for (int i = 0; i < n_samples; ++i) {
      const double tau = t_gen_s + window_s * static_cast<double>(i) / (n_samples - 1);
      Cart6 x = agent->GetStateAt(Real(tau));
      if (i == 0) agent_frame = x.GetFrame();
      t_s(i) = tau;
      rv.row(i) = x.cast<double>().transpose();
    }

    // Broadcast in output_frame: convert the sampled arc from the agent's dynamics
    // frame if needed (tau is the absolute epoch each state was propagated to).
    if (agent_frame != config_.output_frame) {
      for (int i = 0; i < n_samples; ++i) {
        const Vec6 rv_in = rv.row(i).transpose().cast<Real>();
        const Vec6 rv_out = ConvertFrame(Real(t_s(i)), rv_in, agent_frame, config_.output_frame);
        rv.row(i) = rv_out.cast<double>().transpose();
      }
    }

    if (is_ephemeris) {
      GenerateEphemerisFromArc(t_gen_s, t_s, rv);
      ephemeris_msgs_.back().frame = config_.output_frame;
    } else {
      GenerateAlmanacFromArc(t_gen_s, t_s, rv);
      almanac_msgs_.back().frame = config_.output_frame;
    }
  }

  const BroadcastMessage* EphemerisGenApp::Latest(const std::vector<BroadcastMessage>& msgs,
                                                  double t_s) {
    const BroadcastMessage* best = nullptr;
    // Later messages supersede earlier ones, so keep the last window that covers t_s.
    for (const auto& m : msgs) {
      if (m.t_start_s <= t_s && t_s <= m.t_end_s) best = &m;
    }
    return best;
  }

}  // namespace lupnt
