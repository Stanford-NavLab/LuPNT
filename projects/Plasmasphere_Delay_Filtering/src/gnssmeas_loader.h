#pragma once

#include <lupnt/lupnt.h>

#include <Eigen/Dense>
#include <string>
#include <vector>

namespace filtering_sim {
  using namespace lupnt;

  using ConstFreqPrnTuple = std::tuple<std::string, std::string, int>;  // (gnss_const, signal, prn)

  struct GnssMeas {
    int tidx;
    double t_tai;
    double tspan;
    double epoch_t;

    std::string gnss_const;
    std::string signal;

    int sat_id;
    int prn;
    double min_alt;

    double sigma_range;
    double sigma_rangerate;
    double sigma_carrier;
    double cn0;

    Vec6 posvel_tx_ecef;
    Vec6 posvel_tx_mci;
    Vec6 posvel_tx_gcrf;
    Vec6 posvel_tx_pa;

    Vec6 posvel_rx_ecef;
    Vec6 posvel_rx_mci;
    Vec6 posvel_rx_gcrf;
    Vec6 posvel_rx_pa;

    double clockbias_tx;

    Vec6 posvel_tx_ephem_ecef;
    Vec6 posvel_tx_ephem_mci;
    Vec6 posvel_tx_ephem_gcrf;
    Vec6 posvel_tx_ephem_pa;

    double clockbias_ephem;
    // double tecu;
    double tec_delay_m;
    // double second_delay_m;
    // double third_delay_m;
    // double dist_bend_m;
    // double tec_delay_bend_m;
    double total_delay_m;
    double total_carrier_delay_m;
    // double max_sep_line_m;
    // double final_pos_err_m;

    double shapiro_delay;
    double relativistic_delay;
    double ephem_range_error_m;
  };

  class GnssMeasLoader {
  public:
    using Container = std::vector<GnssMeas>;
    using TidxGroups = std::map<int, Container>;

    GnssMeasLoader() = default;

    // First load receiver history
    void load_rx_history(std::filesystem::path cache_path);

    // Load everything in one shot
    void load(const std::string& csv_path, std::filesystem::path cache_path, bool recompute = false,
              bool correct_pco = true, bool correct_clock_bias = true,
              double max_ephem_error_m = 5.0);

    void set_min_cn0_dbhz(double min_cn0_dbhz) { min_cn0_dbhz_ = min_cn0_dbhz; }
    void construct_merged_posvel_histories() { construct_posvel_histories_(data_); }

    void print_delay_statistics();

    // Return only measurements with the given PRN
    static Container filter_by_prn(const Container& all, int prn);
    // Group measurements by tidx (time index)
    static TidxGroups group_by_tidx(const Container& all);

    Container get_data() const { return data_; }
    double get_t0_tai() const { return tais_(0); }
    VecXd get_tspan() const { return tspan_; }
    VecXd get_tais() const { return tais_; }
    MatX6d get_receiver_posvel_ecef() const { return posvel_rx_ecef_; }
    MatX6d get_receiver_posvel_mci() const { return posvel_rx_mci_; }
    std::vector<int> get_dataidx_for_tidx(int tidx) const {
      auto it = tidx_to_dataidx_.find(tidx);
      if (it != tidx_to_dataidx_.end()) {
        return it->second;
      } else {
        return {};
      }
    }
    std::vector<ConstFreqPrnTuple> get_tuples_for_tidx(int tidx) const {
      auto it = tidx_to_const_signal_prn_.find(tidx);
      if (it != tidx_to_const_signal_prn_.end()) {
        return it->second;
      } else {
        return {};
      }
    }

  private:
    double min_cn0_dbhz_ = 18.0;  // [dB-Hz] Minimum C/N0 for measurements

    static std::vector<std::string> split_csv_line(const std::string& line);
    static Vec3 parse_vec3(const std::string& field);
    static std::string trim(const std::string& s);
    Container load_csv_data_(const std::string& csv_path, H5Easy::File& cache_file,
                             bool recompute = false, bool correct_pco = true,
                             bool correct_clock_bias = true, double max_ephem_error_m = 5.0);
    Container convert_frames_(Container& data, H5Easy::File& cache_file);
    Container load_rv_cache_(Container& data, H5Easy::File& cache_file);
    void save_rv_to_cache_(Container& data, H5Easy::File& cache_file);
    void construct_posvel_histories_(Container& data);
    void clear_history_vectors();

    Container load_h5_data_(H5Easy::File& cache_file);
    void save_csv_data_to_cache_(Container& data, H5Easy::File& cache_file);

    // Data
    Container data_;  // loaded container
    VecX tspan_;
    VecX tais_;

    // Position-velocity histories
    int n_steps = 0;
    MatX6 posvel_rx_ecef_;  // receiver position and velocity history
    MatX6 posvel_tx_ecef_;  // transmitter position and velocity history
    MatX6 posvel_tx_mci_;
    MatX6 posvel_rx_mci_;         // receiver position and velocity history in MCI frame
    MatX6 posvel_tx_pa_;          // transmitter position and velocity history in PA frame
    MatX6 posvel_rx_pa_;          // receiver position and velocity history in PA frame
    MatX6 posvel_tx_gcrf_;        // transmitter position and velocity history in GCRF frame
    MatX6 posvel_rx_gcrf_;        // receiver position and velocity history in GCRF frame
    MatX6 posvel_tx_ephem_ecef_;  // transmitter position and velocity history from ephemeris in MCI
                                  // frame
    MatX6 posvel_tx_ephem_mci_;   // transmitter position and velocity history from ephemeris in MCI
                                  // frame
    MatX6 posvel_tx_ephem_gcrf_;  // transmitter position and velocity history from ephemeris in
                                  // ECEF frame
    MatX6 posvel_tx_ephem_pa_;    // transmitter position and velocity history from ephemeris in PA
                                  // frame
    VecX shapiro_delays_;
    VecX relativistic_delays_;

    // Data Matrix
    MatX data_matrix_;

    std::map<int, std::vector<int>> tidx_to_dataidx_;  // mapping from time index to data_ index
    std::map<int, std::vector<ConstFreqPrnTuple>>
        tidx_to_const_signal_prn_;  // mapping from time index to (gnss_const, signal, prn)
  };
}  // namespace filtering_sim
