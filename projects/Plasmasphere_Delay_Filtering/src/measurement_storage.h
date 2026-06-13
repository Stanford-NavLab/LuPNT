#pragma once

#include <lupnt/lupnt.h>

namespace filtering_sim {
  using namespace lupnt;

  enum MeasType {
    CODE = 0,
    CARRIER = 1,
    GRAPHIC = 2,
    TDCP = 3,
    IONOFREE_CODE = 4,
    IONOFREE_CARRIER = 5,
    IONOFREE_TDCP = 6,
  };

  enum MeasFreq {
    L1 = 0,
    L5 = 1,
    L1_L5 = 2,
  };

  enum GnssConst {
    GPS = 0,
    GALILEO = 1,
    QZSS = 2,
    BEIDOU = 3,
    GLONASS = 4,
  };

  enum TrackType {
    ACQ = 0,  // Acquisition
    TRK = 1,  // Tracking
  };

  // Define Freq, PRN pair
  using FreqPrnPair = std::pair<MeasFreq, int>;

  class MeasStorage {
  public:
    std::vector<int> tidx;
    std::vector<double> tspan;
    std::vector<double> tai;
    std::vector<MeasFreq> freq;
    std::vector<int> prn;

    // Measurements
    std::vector<double> pseudorange;
    std::vector<double> carrier_phase_range;
    std::vector<double> carrier_phase;

    // Measurement errors
    std::vector<double> pr_delay;
    std::vector<double> carrier_delay;
    std::vector<int> carrier_phase_integer;
    std::vector<double> sigma_pseudorange;
    std::vector<double> sigma_carrier_phase;

    // Cycle slips
    std::vector<int> cycle_slip;

    // Gnss Parameters
    std::vector<Vec6d> gnss_rv_ephem_mci;
    std::vector<double> gnss_clock_ephem;

    // Tracking status
    std::vector<bool> acquired_this_step;
    std::vector<bool> final_track_this_step;

    // Maps
    std::map<int, std::vector<int>> index_at_tidx;  // indices at each time index
    std::map<int, double> tspan_at_tidx;
    std::map<int, double> tai_at_tidx;
    std::map<int, std::vector<FreqPrnPair>>
        freqs_prns_at_tidx;  // freqs and prns at each time index
    std::map<int, int> num_meas_at_tidx;

    // Unique PRNs
    std::set<int> unique_prns;

    // Packed key -> index
    std::unordered_map<uint64_t, int> index;

    void reserve(size_t n) {
      tidx.reserve(n);
      tspan.reserve(n);
      tai.reserve(n);
      freq.reserve(n);
      prn.reserve(n);
      pseudorange.reserve(n);
      carrier_phase_range.reserve(n);
      carrier_phase.reserve(n);
      carrier_phase_integer.reserve(n);
      pr_delay.reserve(n);
      carrier_delay.reserve(n);
      sigma_pseudorange.reserve(n);
      sigma_carrier_phase.reserve(n);
      gnss_rv_ephem_mci.reserve(n);
      gnss_clock_ephem.reserve(n);
      index.reserve(n);
      cycle_slip.reserve(n);
    }

    void add_measurement(int time_index, double tspan_, double tai_, MeasFreq meas_freq, int prn_,
                         double pr, double carrier_phase_range_, double carrier_phase_,
                         int carrier_phase_integer_ = 0, int cycle_slip_ = 0, double pr_delay = 0.0,
                         double carrier_delay = 0.0, double sigma_pr = 0.0,
                         double sigma_carrier = 0.0, Vec6d gnss_rv_ephem_mci_ = Vec6d::Zero(),
                         double gnss_clock_ephem_ = 0.0, bool acquired_this_step_ = false,
                         bool final_track_this_step_ = false);
    int get_index(int time_index, MeasFreq meas_freq, int prn_) const;
    std::vector<int> get_index_at_tidx(int time_index) const;
    std::vector<std::pair<int, int>> GetDualFreqPairs(int time_index) const;
    int GetUniquePrnCount() const { return static_cast<int>(unique_prns.size()); }
    std::vector<int> GetNewAcquiredSignalsAtTimeIndex(int time_index) const;
    std::vector<int> GetFinalTrackSignalsAtTimeIndex(int time_index) const;

  private:
    // [ 32 bits tidx ][ 10 bits prn ][ 2 bits freq ]
    static inline uint64_t make_key(int time_index, MeasFreq f, int prn_);
  };

  // Other functions related to measurement storage can be declared here
  int GetUniquePrn(int prn, GnssConst gnss_const);
  ClockModel StringToClockModel(const std::string& model_str);
  GnssConst StringToGnssConst(const std::string& str);
  MeasFreq StringToMeasFreq(const std::string& str);
  double get_wavelength(MeasFreq freq);

}  // namespace filtering_sim
