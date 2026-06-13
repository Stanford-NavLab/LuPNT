

#include "src/measurement_storage.h"

#include <lupnt/lupnt.h>

#include <cassert>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace filtering_sim {
  using namespace lupnt;
  template <typename T> using vec = std::vector<T>;

  /***************************************************************************
   * Measurement Storage Class
   * **************************************************************************/
  void MeasStorage::add_measurement(int time_index, double tspan_, double tai_, MeasFreq meas_freq,
                                    int prn_, double pr, double carrier_phase_range_,
                                    double carrier_phase_, int carrier_phase_integer_,
                                    int cycle_slip_, double pr_delay_, double carrier_delay_,
                                    double sigma_pr, double sigma_carrier, Vec6d gnss_rv_ephem_mci_,
                                    double gnss_clock_ephem_, bool acquired_this_step_,
                                    bool final_track_this_step_) {
    const int i = static_cast<int>(tidx.size());

    tidx.push_back(time_index);
    tspan.push_back(tspan_);
    tai.push_back(tai_);
    freq.push_back(meas_freq);
    prn.push_back(prn_);
    pseudorange.push_back(pr);
    carrier_phase_range.push_back(carrier_phase_range_);
    carrier_phase.push_back(carrier_phase_);
    carrier_phase_integer.push_back(carrier_phase_integer_);
    cycle_slip.push_back(cycle_slip_);
    sigma_pseudorange.push_back(sigma_pr);
    sigma_carrier_phase.push_back(sigma_carrier);
    gnss_rv_ephem_mci.push_back(gnss_rv_ephem_mci_);
    gnss_clock_ephem.push_back(gnss_clock_ephem_);
    pr_delay.push_back(pr_delay_);
    carrier_delay.push_back(carrier_delay_);
    acquired_this_step.push_back(acquired_this_step_);
    final_track_this_step.push_back(final_track_this_step_);

    // Update index_at_tidx map
    index_at_tidx[time_index].push_back(i);
    num_meas_at_tidx[time_index] += 1;

    // PRNs at time index
    freqs_prns_at_tidx[time_index].push_back({meas_freq, prn_});
    // tspan and tai at time index
    tspan_at_tidx[time_index] = tspan_;
    tai_at_tidx[time_index] = tai_;

    // Update unique PRNs set
    unique_prns.insert(prn_);

    // Last-write-wins for duplicate keys
    index[make_key(time_index, meas_freq, prn_)] = i;
  }

  int MeasStorage::get_index(int time_index, MeasFreq meas_freq, int prn_) const {
    const uint64_t k = make_key(time_index, meas_freq, prn_);
    auto it = index.find(k);
    return (it == index.end()) ? -1 : it->second;
  }

  std::vector<int> MeasStorage::get_index_at_tidx(int time_index) const {
    auto it = index_at_tidx.find(time_index);
    if (it == index_at_tidx.end()) {
      return {};
    } else {
      return it->second;
    }
  }

  uint64_t MeasStorage::make_key(int time_index, MeasFreq f, int prn_) {
    assert(prn_ >= 0 && prn_ <= 1023);
    const uint64_t ti = static_cast<uint32_t>(time_index);
    const uint64_t pr = static_cast<uint32_t>(prn_);
    const uint64_t fr = static_cast<uint8_t>(f);  // 0..2

    return (ti << 12) | (pr << 2) | fr;
  }

  std::vector<std::pair<int, int>> MeasStorage::GetDualFreqPairs(int time_index) const {
    // Return the vector of (L1_index, L5_index) pairs for dual-frequency measurements at the given
    // time index The pair shares prn but differs in frequency (L1, L5)
    std::vector<std::pair<int, int>> dual_freq_pairs;

    // First, build a map from prn to indices for L1 and L5
    std::map<int, std::pair<int, int>> prn_to_idx_L1L5;
    std::vector<int> indices = get_index_at_tidx(time_index);
    for (int idx : indices) {
      MeasFreq f = freq[idx];
      int prn_ = prn[idx];

      // Initialize if not present
      if (prn_to_idx_L1L5.find(prn_) == prn_to_idx_L1L5.end()) {
        prn_to_idx_L1L5[prn_] = {-1, -1};
      }

      if (f == MeasFreq::L1) {
        prn_to_idx_L1L5[prn_].first = idx;  // L1 index
      } else if (f == MeasFreq::L5) {
        prn_to_idx_L1L5[prn_].second = idx;  // L5 index
      }
    }

    // Collect pairs where both L1 and L5 measurements exist
    for (const auto& kv : prn_to_idx_L1L5) {
      if (kv.second.first != -1 && kv.second.second != -1) {
        dual_freq_pairs.push_back(kv.second);
      }
    }

    return dual_freq_pairs;
  }

  std::vector<int> MeasStorage::GetNewAcquiredSignalsAtTimeIndex(int time_index) const {
    std::vector<int> acquired_indices;  // Indices of newly acquired signals at the given time index
    std::vector<FreqPrnPair> freq_prn_pairs = freqs_prns_at_tidx.at(time_index);
    std::vector<FreqPrnPair> prev_freq_prn_pairs = freqs_prns_at_tidx.at(time_index - 1);
    std::set<FreqPrnPair> prev_set(prev_freq_prn_pairs.begin(), prev_freq_prn_pairs.end());
    for (const auto& fp : freq_prn_pairs) {
      if (prev_set.find(fp) == prev_set.end()) {
        int idx = get_index(time_index, fp.first, fp.second);
        if (idx != -1) {
          acquired_indices.push_back(idx);
        }
      }
    }

    return acquired_indices;
  }

  std::vector<int> MeasStorage::GetFinalTrackSignalsAtTimeIndex(int time_index) const {
    std::vector<int> final_track_indices;  // Indices of signals tracked at the final time index
    std::vector<FreqPrnPair> freq_prn_pairs = freqs_prns_at_tidx.at(time_index);
    std::vector<FreqPrnPair> next_freq_prn_pairs = freqs_prns_at_tidx.at(time_index + 1);
    std::set<FreqPrnPair> next_set(next_freq_prn_pairs.begin(), next_freq_prn_pairs.end());
    for (const auto& fp : freq_prn_pairs) {
      if (next_set.find(fp) == next_set.end()) {
        // Not found in next time index => final track
        int idx = get_index(time_index, fp.first, fp.second);
        if (idx != -1) {
          final_track_indices.push_back(idx);
        }
      }
    }

    return final_track_indices;
  }

  /***************************************************************************
   * Additional Measurement Storage Utilities
   * **************************************************************************/

  int GetUniquePrn(int prn, GnssConst gnss_const) {
    int offset = 0;
    switch (gnss_const) {
      case GPS: offset = 0; break;
      case GALILEO: offset = 40; break;
      case QZSS: offset = 80; break;
      case BEIDOU: offset = 100; break;
      case GLONASS: offset = 150; break;
      default: offset = 0; break;
    }
    return prn + offset;
  }

  ClockModel StringToClockModel(const std::string& model_str) {
    if (model_str == "OCXO") {
      return ClockModel::OCXO;
    } else if (model_str == "USO") {
      return ClockModel::USO;
    } else if (model_str == "CSAC") {
      return ClockModel::CSAC;
    } else if (model_str == "MINI_RAFS") {
      return ClockModel::MINI_RAFS;
    } else if (model_str == "RAFS") {
      return ClockModel::RAFS;
    } else if (model_str == "DSAC") {
      return ClockModel::DSAC;
    } else {
      return ClockModel::UNDEFINED;
    }
  }

  GnssConst StringToGnssConst(const std::string& str) {
    if ((str == "GPS") || (str == "Gps") || (str == "gps")) {
      return GnssConst::GPS;
    } else if ((str == "GALILEO") || (str == "Galileo") || (str == "galileo")) {
      return GnssConst::GALILEO;
    } else if ((str == "QZSS") || (str == "Qzss") || (str == "qzss")) {
      return GnssConst::QZSS;
    } else if ((str == "BEIDOU") || (str == "Beidou") || (str == "beidou")) {
      return GnssConst::BEIDOU;
    } else if ((str == "GLONASS") || (str == "Glonass") || (str == "glonass")) {
      return GnssConst::GLONASS;
    } else {
      throw std::invalid_argument("Unknown GNSS constellation string: " + str);
    }
  }

  MeasFreq StringToMeasFreq(const std::string& str) {
    if ((str == "L1") || (str == "E1")) {
      return MeasFreq::L1;
    } else if ((str == "L5") || (str == "E5a")) {
      return MeasFreq::L5;
    } else if ((str == "L1_L5") || (str == "E1_E5a")) {
      return MeasFreq::L1_L5;
    } else {
      throw std::invalid_argument("Unknown measurement frequency string: " + str);
    }
  }

  double get_wavelength(MeasFreq freq) {
    double FREQ_L1 = 1575.42e6;  // [Hz]
    double FREQ_L5 = 1176.45e6;  // [Hz]
    double FREQ_L1_L5 = 0.0;     // Not a physical frequency

    double iono_factor_L1 = 2.2606043275188257;
    double iono_factor_L5 = 1.2606043275188257;

    double lambda_L1 = 0.19029367279836487;  // C / FREQ_L1
    double lambda_L5 = 0.25482804879085386;  // C / FREQ_L5
    double lambda_L1_L5
        = iono_factor_L1 * lambda_L1 - iono_factor_L5 * lambda_L5;  // Not a physical wavelength

    switch (freq) {
      case MeasFreq::L1: return lambda_L1;
      case MeasFreq::L5: return lambda_L5;
      case MeasFreq::L1_L5: return lambda_L1_L5;
      default: throw std::invalid_argument("Unsupported frequency for wavelength calculation.");
    }
  }

}  // namespace filtering_sim
