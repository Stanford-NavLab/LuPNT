#pragma once

#include <utility>
#include <vector>

#include "lupnt/core/definitions.h"
#include "lupnt/measurements/gnss_measurement.h"
#include "lupnt/measurements/measurement.h"
#include "lupnt/measurements/measurement_utils.h"

namespace lupnt {

  /// @brief Combined lunar-GNSS ODTS measurement model: a stack of per-channel GNSS
  /// observables (pseudorange/Doppler/carrier-phase) followed, when time-differenced carrier
  /// phase (TDCP) is enabled, by TDCP rows built from a stochastic-cloning state
  /// `[x_current; x_previous]`.
  ///
  /// This owns the measurement-model math that previously lived inside
  /// `LunarGnssODTSSimulation`: the per-channel loop delegates to `GnssMeasurement`, and the
  /// TDCP rows difference the carrier range between the current and previous cloned states.
  /// The reusable pieces are exposed as static helpers so the simulation's truth-data
  /// generation shares exactly the same model.
  class LunarGnssCombinedMeasurement : public MeasurementClone<LunarGnssCombinedMeasurement> {
  public:
    /// @brief A matched carrier-phase channel pair (same constellation/PRN/frequency) across
    /// two consecutive epochs, used to form a TDCP observable.
    struct TdcpPair {
      GnssChannel current;
      GnssChannel previous;
    };

    struct Config {
      std::vector<GnssChannel> channels;           ///< current-epoch channels (pseudorange/Doppler)
      GnssMeasurementOptions current_options;      ///< observable/state config for `channels`
      GnssMeasurementOptions carrier_options;      ///< carrier-phase config for TDCP rows
      std::vector<TdcpPair> tdcp_pairs;            ///< matched carrier-phase pairs (TDCP)
      bool use_tdcp = false;                       ///< enable TDCP rows / cloned-state layout
      double tdcp_sigma_m = 0.0;                   ///< shared additive TDCP sigma [m]
      double filter_tdcp_noise_inflation_m = 0.0;  ///< filter-only TDCP inflation [m]
    };

    LunarGnssCombinedMeasurement() = default;
    explicit LunarGnssCombinedMeasurement(Config config) : config_(std::move(config)) {}

    const Config& GetConfig() const { return config_; }
    void SetConfig(const Config& config) { config_ = config; }

    MeasData Compute(const State& x, MatXd* H = nullptr) const override;

    // ---- Reusable model helpers (shared with truth-data generation) ----------------------

    /// @brief Stack the per-channel GNSS observable vector `y = h(x)` (and, if `H != nullptr`,
    /// its Jacobian) by delegating to `GnssMeasurement` for each channel.
    static VecXd ComputeMeasurementVector(const State& x, const std::vector<GnssChannel>& channels,
                                          const GnssMeasurementOptions& options,
                                          MatXd* H = nullptr);

    /// @brief Diagonal measurement covariance built from each channel's per-observable sigmas.
    static MatXd MeasurementCovariance(const std::vector<GnssChannel>& channels,
                                       const GnssMeasurementOptions& options);

    /// @brief Match carrier-phase channels between two epochs (same const/PRN/frequency).
    static std::vector<TdcpPair> MakeTdcpPairs(const std::vector<GnssChannel>& current,
                                               const std::vector<GnssChannel>& previous);

    /// @brief TDCP measurement vector (per pair, current-minus-previous carrier range) and,
    /// if `H != nullptr`, its Jacobian with respect to `[x_current; x_previous]`.
    static VecXd ComputeTdcpVector(const State& current_state, const State& previous_state,
                                   const std::vector<TdcpPair>& pairs,
                                   const GnssMeasurementOptions& options, MatXd* H = nullptr);

    /// @brief Diagonal TDCP covariance. `filter_tdcp_inflation_m` is added in quadrature
    /// (pass 0 for the truth-side covariance, the filter inflation for the filter side).
    static MatXd TdcpCovariance(const std::vector<TdcpPair>& pairs, double tdcp_sigma_m,
                                double filter_tdcp_inflation_m);

  private:
    /// @brief Carrier-phase range in meters for one channel (and its Jacobian if requested).
    static Real CarrierRangeMeters(const State& x, const GnssChannel& channel,
                                   const GnssMeasurementOptions& options, MatXd* H = nullptr);

    Config config_;
  };

}  // namespace lupnt
