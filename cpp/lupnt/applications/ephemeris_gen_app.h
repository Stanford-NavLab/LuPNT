#pragma once

#include <vector>

#include "lupnt/applications/lunanet_almanac.h"
#include "lupnt/applications/lunanet_ephemeris.h"
#include "lupnt/applications/lunanet_sat_app.h"
#include "lupnt/conversions/frame_converter.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief One generated broadcast navigation message (ephemeris or almanac),
  /// valid over a single fitting window.
  ///
  /// `params` is the fitted parameter vector of the corresponding model
  /// (`CartesianEphemeris`/`Almanac`); see that model's `ParamNames()` for the
  /// layout. In a real routine these are the bits that would be downlinked and
  /// re-evaluated by a receiver over `[t_start_s, t_end_s]`.
  struct BroadcastMessage {
    double t_generated_s = 0.0;    // elapsed time this message was generated [s]
    double t_start_s = 0.0;        // validity-window start [s, same origin]
    double t_end_s = 0.0;          // validity-window end [s, same origin]
    Frame frame = Frame::MOON_CI;  // frame of the fitted arc (= agent dynamics frame)
    VecXd params;                  // fitted parameter vector
  };

  /// @brief Configuration for `EphemerisGenApp`.
  struct EphemerisGenConfig {
    /// Generate a precise, short-validity broadcast ephemeris (`CartesianEphemeris`).
    bool generate_ephemeris = true;
    /// Generate a coarse, long-validity broadcast almanac (`Almanac`).
    bool generate_almanac = true;

    /// Short-validity ephemeris fit options and its validity window / refresh
    /// cadence. `ephemeris_refresh_s <= 0` refreshes once per validity window.
    EphemerisFitOptions ephemeris_options;
    double ephemeris_window_s = 2.0 * SECS_HOUR;
    double ephemeris_refresh_s = 0.0;

    /// Long-validity almanac fit options and its validity window / refresh
    /// cadence. `almanac_refresh_s <= 0` refreshes once per validity window.
    AlmanacFitOptions almanac_options;
    double almanac_window_s = 15.0 * SECS_DAY;
    double almanac_refresh_s = 0.0;

    /// Number of samples of the predicted arc used to fit each window (queried
    /// uniformly across the validity window via the owning agent's
    /// `GetStateAt`). Should comfortably exceed the polynomial order.
    int ephemeris_fit_samples = 121;
    int almanac_fit_samples = 361;

    /// Frame the broadcast messages are fit and evaluated in. When this differs
    /// from the agent's dynamics frame, the predicted arc sampled from the agent
    /// is converted to `output_frame` before fitting (e.g. `Frame::MOON_PA` to
    /// broadcast in the rotating Moon-fixed principal-axis frame, per Iiyama &
    /// Gao). `Setup` forces both `ephemeris_options.frame` and
    /// `almanac_options.frame` to this value.
    Frame output_frame = Frame::MOON_CI;
  };

  /// @brief LunaNet sub-app that generates broadcast ephemeris (precise,
  /// short-validity) and almanac (coarse, long-validity) navigation messages
  /// from the satellite's own predicted trajectory.
  ///
  /// This is the "navigation-message generation" sub-app anticipated by
  /// `LunaNetSubApp`: attach it to a `LunaNetSatApp` via `AddSubApp` to embed
  /// broadcast-message generation in the LunaNet satellite routine. On each
  /// `Step(t)`, when a refresh is due, it samples the owning agent's predicted
  /// state over the corresponding validity window (`Agent::GetStateAt`, i.e. it
  /// propagates its own orbit forward the way a real satellite would before
  /// broadcasting) and fits a `CartesianEphemeris` / `Almanac`, appending the
  /// resulting `BroadcastMessage`.
  ///
  /// When no agent is attached (e.g. unit tests or offline generation) the same
  /// models can be driven directly via `GenerateEphemerisFromArc` /
  /// `GenerateAlmanacFromArc`.
  class EphemerisGenApp : public LunaNetSubApp {
  public:
    explicit EphemerisGenApp(EphemerisGenConfig config = {});

    void Setup(LunaNetSatApp& app) override;
    void Step(Real t) override;
    void Finish() override;

    /// @brief Fit and store an ephemeris message from an explicit sampled arc.
    /// @param t_gen_s  Elapsed time the message is generated / becomes valid [s].
    /// @param t_s      Sample epochs of the arc [s] (strictly increasing).
    /// @param rv       Sampled Cartesian states [N x 6].
    const BroadcastMessage& GenerateEphemerisFromArc(double t_gen_s, const VecXd& t_s,
                                                     const MatXd& rv);
    /// @brief Fit and store an almanac message from an explicit sampled arc.
    const BroadcastMessage& GenerateAlmanacFromArc(double t_gen_s, const VecXd& t_s,
                                                   const MatXd& rv);

    const std::vector<BroadcastMessage>& GetEphemerisMessages() const { return ephemeris_msgs_; }
    const std::vector<BroadcastMessage>& GetAlmanacMessages() const { return almanac_msgs_; }

    /// @brief Most recent generated message whose validity window contains
    /// `t_s` (nullptr if none) -- mirrors a receiver selecting the currently
    /// valid broadcast page.
    const BroadcastMessage* LatestEphemeris(double t_s) const;
    const BroadcastMessage* LatestAlmanac(double t_s) const;

    const EphemerisGenConfig& GetConfig() const { return config_; }
    const CartesianEphemeris& GetEphemerisModel() const { return ephemeris_; }
    const Almanac& GetAlmanacModel() const { return almanac_; }

  private:
    // Sample the owning agent's predicted arc over [t_gen, t_gen + window] at
    // `n_samples` uniform epochs and fit the ephemeris (`is_ephemeris`) or almanac.
    void GenerateFromAgent(double t_gen_s, double window_s, int n_samples, bool is_ephemeris);
    static const BroadcastMessage* Latest(const std::vector<BroadcastMessage>& msgs, double t_s);

    EphemerisGenConfig config_;
    CartesianEphemeris ephemeris_;
    Almanac almanac_;
    LunaNetSatApp* sat_app_ = nullptr;

    std::vector<BroadcastMessage> ephemeris_msgs_;
    std::vector<BroadcastMessage> almanac_msgs_;
    double next_ephemeris_gen_s_ = 0.0;
    double next_almanac_gen_s_ = 0.0;
  };

}  // namespace lupnt
