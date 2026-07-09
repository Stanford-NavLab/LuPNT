#include "lupnt/measurements/lunar_gnss_combined_measurement.h"

#include <algorithm>
#include <cmath>

namespace lupnt {

  VecXd LunarGnssCombinedMeasurement::ComputeMeasurementVector(
      const State& x, const std::vector<GnssChannel>& channels,
      const GnssMeasurementOptions& options, MatXd* H) {
    const int n_obs = static_cast<int>(options.observables.size());
    VecXd y(n_obs * static_cast<int>(channels.size()));
    if (H != nullptr) H->setZero(y.size(), x.size());

    for (int i = 0; i < static_cast<int>(channels.size()); ++i) {
      GnssMeasurement measurement(channels[i]);
      MatXd H_i;
      VecXd y_i = measurement.ComputeVector(x, H != nullptr ? &H_i : nullptr, options);
      y.segment(i * n_obs, n_obs) = y_i;
      if (H != nullptr) H->block(i * n_obs, 0, n_obs, x.size()) = H_i;
    }
    return y;
  }

  MatXd LunarGnssCombinedMeasurement::MeasurementCovariance(
      const std::vector<GnssChannel>& channels, const GnssMeasurementOptions& options) {
    const int n_obs = static_cast<int>(options.observables.size());
    MatXd R = MatXd::Zero(n_obs * static_cast<int>(channels.size()),
                          n_obs * static_cast<int>(channels.size()));
    for (int i = 0; i < static_cast<int>(channels.size()); ++i) {
      for (int j = 0; j < n_obs; ++j) {
        const int row = i * n_obs + j;
        if (options.observables[j] == GnssObservable::PSEUDORANGE)
          R(row, row) = std::pow(channels[i].sigma_pseudorange_m.val(), 2);
        else if (options.observables[j] == GnssObservable::DOPPLER)
          R(row, row) = std::pow(channels[i].sigma_doppler_hz.val(), 2);
        else if (options.observables[j] == GnssObservable::CARRIER_PHASE)
          R(row, row) = std::pow(channels[i].sigma_carrier_phase_cycles.val(), 2);
      }
    }
    return R;
  }

  std::vector<LunarGnssCombinedMeasurement::TdcpPair> LunarGnssCombinedMeasurement::MakeTdcpPairs(
      const std::vector<GnssChannel>& current, const std::vector<GnssChannel>& previous) {
    std::vector<TdcpPair> pairs;
    for (const auto& curr : current) {
      auto it = std::find_if(previous.begin(), previous.end(), [&](const GnssChannel& prev) {
        return curr.gnss_const == prev.gnss_const && curr.prn == prev.prn
               && curr.frequency == prev.frequency;
      });
      if (it != previous.end()) pairs.push_back({curr, *it});
    }
    return pairs;
  }

  Real LunarGnssCombinedMeasurement::CarrierRangeMeters(const State& x, const GnssChannel& channel,
                                                        const GnssMeasurementOptions& options,
                                                        MatXd* H) {
    GnssMeasurement measurement(channel);
    MatXd H_cycles;
    VecXd y_cycles = measurement.ComputeVector(x, H != nullptr ? &H_cycles : nullptr, options);
    const Real lambda = channel.Wavelength();
    if (H != nullptr) *H = lambda.val() * H_cycles;
    return lambda * y_cycles(0);
  }

  VecXd LunarGnssCombinedMeasurement::ComputeTdcpVector(const State& current_state,
                                                        const State& previous_state,
                                                        const std::vector<TdcpPair>& pairs,
                                                        const GnssMeasurementOptions& options,
                                                        MatXd* H) {
    VecXd y(pairs.size());
    if (H != nullptr) H->setZero(pairs.size(), current_state.size() + previous_state.size());
    for (int i = 0; i < static_cast<int>(pairs.size()); ++i) {
      MatXd H_curr;
      MatXd H_prev;
      Real curr_m = CarrierRangeMeters(current_state, pairs[i].current, options,
                                       H != nullptr ? &H_curr : nullptr);
      Real prev_m = CarrierRangeMeters(previous_state, pairs[i].previous, options,
                                       H != nullptr ? &H_prev : nullptr);
      y(i) = (curr_m - prev_m).val();
      if (H != nullptr) {
        H->block(i, 0, 1, current_state.size()) = H_curr;
        H->block(i, current_state.size(), 1, previous_state.size()) = -H_prev;
      }
    }
    return y;
  }

  MatXd LunarGnssCombinedMeasurement::TdcpCovariance(const std::vector<TdcpPair>& pairs,
                                                     double tdcp_sigma_m,
                                                     double filter_tdcp_inflation_m) {
    const double filter_inflation2 = std::pow(filter_tdcp_inflation_m, 2);
    MatXd R = MatXd::Zero(pairs.size(), pairs.size());
    for (int i = 0; i < static_cast<int>(pairs.size()); ++i) {
      const double curr_sigma_m
          = pairs[i].current.sigma_carrier_phase_cycles.val() * pairs[i].current.Wavelength().val();
      const double prev_sigma_m = pairs[i].previous.sigma_carrier_phase_cycles.val()
                                  * pairs[i].previous.Wavelength().val();
      R(i, i) = std::pow(curr_sigma_m, 2) + std::pow(prev_sigma_m, 2) + std::pow(tdcp_sigma_m, 2)
                + filter_inflation2;
    }
    return R;
  }

  MeasData LunarGnssCombinedMeasurement::Compute(const State& x, MatXd* H) const {
    const Config& cfg = config_;
    const int base_n = cfg.use_tdcp ? static_cast<int>(x.size() / 2) : static_cast<int>(x.size());
    State x_curr = cfg.use_tdcp ? State(x.head(base_n)) : x;
    State x_prev = cfg.use_tdcp ? State(x.tail(base_n)) : x;

    MatXd H_current_base;
    VecXd y_current = ComputeMeasurementVector(x_curr, cfg.channels, cfg.current_options,
                                               H != nullptr ? &H_current_base : nullptr);
    MatXd H_tdcp;
    VecXd y_tdcp = cfg.use_tdcp
                       ? ComputeTdcpVector(x_curr, x_prev, cfg.tdcp_pairs, cfg.carrier_options,
                                           H != nullptr ? &H_tdcp : nullptr)
                       : VecXd(0);

    VecXd y(y_current.size() + y_tdcp.size());
    y << y_current, y_tdcp;

    if (H != nullptr) {
      H->setZero(y.size(), x.size());
      H->block(0, 0, y_current.size(), base_n) = H_current_base;
      if (cfg.use_tdcp && y_tdcp.size() > 0) {
        H->block(y_current.size(), 0, y_tdcp.size(), 2 * base_n) = H_tdcp;
      }
    }

    MeasData md;
    md.value = y;
    MatXd R_current = MeasurementCovariance(cfg.channels, cfg.current_options);
    MatXd R_tdcp = cfg.use_tdcp ? TdcpCovariance(cfg.tdcp_pairs, cfg.tdcp_sigma_m,
                                                 cfg.filter_tdcp_noise_inflation_m)
                                : MatXd::Zero(0, 0);
    md.covariance = MatXd::Zero(y.size(), y.size());
    if (R_current.size() > 0)
      md.covariance.topLeftCorner(R_current.rows(), R_current.cols()) = R_current;
    if (R_tdcp.size() > 0) md.covariance.bottomRightCorner(R_tdcp.rows(), R_tdcp.cols()) = R_tdcp;
    return md;
  }

}  // namespace lupnt
