#include "gnssmeas_loader.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "delays.h"
#include "highfive/H5Easy.hpp"
#include "lupnt/lupnt.h"

namespace {
  // simple trim helpers
  inline void ltrim_inplace(std::string& s) {
    s.erase(s.begin(),
            std::find_if(s.begin(), s.end(), [](unsigned char ch) { return !std::isspace(ch); }));
  }

  inline void rtrim_inplace(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); })
                .base(),
            s.end());
  }

}  // namespace

namespace filtering_sim {

  using namespace lupnt;

  // Converters
  double gnss_const_to_double(const std::string& gnss_const) {
    if (gnss_const == "GPS") return 0;
    if (gnss_const == "GLONASS") return 1;
    if (gnss_const == "GALILEO") return 2;
    if (gnss_const == "BEIDOU") return 3;
    if (gnss_const == "QZSS") return 4;
    throw std::runtime_error("Unknown GNSS constellation: " + gnss_const);
  }

  double signal_to_double(const std::string& signal) {
    if (signal == "L1" || signal == "E1") return 0;
    if (signal == "L2") return 1;
    if (signal == "L5" || signal == "E5a") return 2;
    if (signal == "E5b") return 3;
    throw std::runtime_error("Unknown GNSS signal: " + signal);
  }

  std::string gnss_const_to_string(double gnss_const) {
    switch (static_cast<int>(gnss_const)) {
      case 0: return "GPS";
      case 1: return "GLONASS";
      case 2: return "GALILEO";
      case 3: return "BEIDOU";
      case 4: return "QZSS";
      default:
        throw std::runtime_error("Unknown GNSS constellation integer: "
                                 + std::to_string(gnss_const));
    }
  }

  std::string signal_to_string(double signal) {
    switch (static_cast<int>(signal)) {
      case 0: return "L1";
      case 1: return "L2";
      case 2: return "L5";
      case 3: return "E5b";
      default: throw std::runtime_error("Unknown GNSS signal integer: " + std::to_string(signal));
    }
  }

  std::string GnssMeasLoader::trim(const std::string& s) {
    std::string out = s;
    ltrim_inplace(out);
    rtrim_inplace(out);
    return out;
  }

  // Simple CSV splitter assuming no commas inside fields
  // (your bracketed vectors use spaces, so this is safe)
  std::vector<std::string> GnssMeasLoader::split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string item;

    while (std::getline(ss, item, ',')) {
      fields.push_back(item);
    }
    return fields;
  }

  Vec3 GnssMeasLoader::parse_vec3(const std::string& field) {
    std::string s = trim(field);

    if (s.empty()) throw std::runtime_error("Empty vector field");

    if (s.front() != '[' || s.back() != ']') {
      throw std::runtime_error("Vector field not in [x y z] format: " + s);
    }

    // strip [ and ]
    s = s.substr(1, s.size() - 2);
    std::stringstream ss(s);

    double x, y, z;
    if (!(ss >> x >> y >> z)) {
      throw std::runtime_error("Failed to parse 3 components from: " + s);
    }
    return Vec3{x, y, z};
  }

  void GnssMeasLoader::load_rx_history(std::filesystem::path cache_path) {
    bool recompute = false;

    H5Easy::File rx_cache_file(cache_path.string(), H5Easy::File::ReadOnly);

    std::cout << "Loading receiver history from cache..." << std::endl;

    // Dummy lambda for compatibility
    auto func = []() { return MatX6(); };
    auto func_vec = []() { return VecX(); };

    tspan_ = LoadOrRecompute<-1, 1, double>("/tspan", rx_cache_file, recompute, func_vec);
    tais_ = LoadOrRecompute<-1, 1, double>("/t_tai", rx_cache_file, recompute, func_vec);
    posvel_rx_ecef_
        = LoadOrRecompute<-1, 6, double>("/posvel_rx_ecef", rx_cache_file, recompute, func);
    posvel_rx_mci_
        = LoadOrRecompute<-1, 6, double>("/posvel_rx_mci", rx_cache_file, recompute, func);
    posvel_rx_gcrf_
        = LoadOrRecompute<-1, 6, double>("/posvel_rx_gcrf", rx_cache_file, recompute, func);
    posvel_rx_pa_ = LoadOrRecompute<-1, 6, double>("/posvel_rx_pa", rx_cache_file, recompute, func);

    n_steps = static_cast<int>(tspan_.size());
  }

  GnssMeasLoader::Container GnssMeasLoader::load_csv_data_(const std::string& csv_path,
                                                           H5Easy::File& cache_file, bool recompute,
                                                           bool correct_pco,
                                                           bool correct_clock_bias,
                                                           double max_ephem_error_m) {
    std::ifstream in(csv_path);

    if (!in.is_open()) {
      throw std::runtime_error("Failed to open GNSS measurement CSV: " + csv_path);
    }

    if (n_steps == 0) {
      throw std::runtime_error("Receiver history not loaded before GNSS measurement loading.");
    }

    std::string line;

    // skip header
    if (!std::getline(in, line)) {
      throw std::runtime_error("GNSS CSV seems empty: " + csv_path);
    }

    // For temp storage of the parameters
    Container data;
    std::cout << "  csv path: " << csv_path << std::endl;

    int rx_history_idx = 0;

    int data_idx = 0;

    std::cout << "  Loading CSV..." << std::endl;

    while (std::getline(in, line)) {
      if (line.empty()) continue;

      auto fields = split_csv_line(line);

      // Expecting exactly 34 columns
      if (fields.size() != 34) {
        std::ostringstream oss;
        oss << "Expected 34 fields, got " << fields.size() << " in line: " << line;
        // skip malformed line
        continue;
      }

      GnssMeas m;
      // For readability, use a small lambda to pull cells
      auto as_int = [&](size_t idx) {
        std::string s = trim(fields[idx]);
        return std::stoi(s);
      };
      auto as_double = [&](size_t idx) {
        std::string s = trim(fields[idx]);
        if (s.empty()) {
          return 0.0;
        }
        return std::stod(s);
      };
      auto as_string = [&](size_t idx) { return trim(fields[idx]); };

      // constants
      const double m_km = 1000.0;       // km to m
      const double c_ms = 299792458.0;  // speed of light in m/s

      // Data
      // tidx,t_tai,tspan,epoch_t,gnss_const,signal,sat_id,prn,min_alt,
      // sigma_range,sigma_rangerate,sigma_carrier,cn0,
      // pos_tx,pos_rx,vel_tx,vel_rx,clockbias_tx,pos_tx_ephem,vel_tx_ephem,
      // clockbias_ephem,tecu,tec_delay_m,second_delay_m,third_delay_m,dist_bend_m,
      // tec_delay_bend_m,total_delay_m,max_sep_line_m,final_pos_err_m
      m.tidx = as_int(0);                                                // 0
      m.t_tai = as_double(1);                                            // 1
      m.tspan = as_double(2);                                            // 2
      m.epoch_t = as_double(3);                                          // 3
      m.gnss_const = as_string(4);                                       // 4
      m.signal = as_string(5);                                           // 5
      m.sat_id = as_int(6);                                              // 6
      m.prn = as_int(7);                                                 // 7
      m.min_alt = as_double(8);                                          // 8
      m.sigma_range = as_double(9) * m_km;                               // 9
      m.sigma_rangerate = as_double(10) * m_km;                          // 10
      m.sigma_carrier = as_double(11) * m_km;                            // 11
      m.cn0 = as_double(12);                                             // 12
      m.posvel_tx_ecef.head<3>() = parse_vec3(fields[13]) * m_km;        // 13
      m.posvel_rx_ecef.head<3>() = parse_vec3(fields[14]) * m_km;        // 14
      m.posvel_tx_ecef.tail<3>() = parse_vec3(fields[15]) * m_km;        // 15
      m.posvel_rx_ecef.tail<3>() = parse_vec3(fields[16]) * m_km;        // 16
      m.clockbias_tx = as_double(17) * C;                                // 17
      m.posvel_tx_ephem_ecef.head<3>() = parse_vec3(fields[18]) * m_km;  // 18
      m.posvel_tx_ephem_ecef.tail<3>() = parse_vec3(fields[19]) * m_km;  // 19
      m.clockbias_ephem = as_double(20) * C;                             // 20
      // m.tecu = as_double(21);                                            // 21
      m.tec_delay_m = as_double(22);          // 22
      double second_delay_m = as_double(23);  // 23
      double third_delay_m = as_double(24);   // 24
      double dist_bend_m = as_double(25);     // 25
      // m.tec_delay_bend_m = as_double(26);                                // 26
      // m.total_delay_m = as_double(27);                                   // 27
      // m.max_sep_line_m = as_double(28);                                  // 28
      // m.final_pos_err_m = as_double(29);                                 // 29
      double clock_bias_tx_median = as_double(30) * C;  // 30
      double pco_ecef_x_m = as_double(31);              // 31
      double pco_ecef_y_m = as_double(32);              // 32
      double pco_ecef_z_m = as_double(33);              // 33

      // Modify the true posvel to include the PCO offset
      if (correct_pco) {
        // Change true side
        // m.posvel_tx_ecef.head<3>() += Vec3{pco_ecef_x_m, pco_ecef_y_m, pco_ecef_z_m};
        // Change ephemeris side
        m.posvel_tx_ephem_ecef.head<3>() -= Vec3{pco_ecef_x_m, pco_ecef_y_m, pco_ecef_z_m};
      }
      if (correct_clock_bias) {
        // Change true side
        // m.clockbias_tx += clock_bias_tx_median;
        // Change ephemeris side
        m.clockbias_ephem -= clock_bias_tx_median;
      }

      if (m.cn0 < min_cn0_dbhz_) {
        continue;  // Ignore low C/N0 measurements
      }

      // Additional measurement: total pseudorange delay and carrier delay
      m.total_delay_m = dist_bend_m + m.tec_delay_m + second_delay_m + third_delay_m;

      if (m.total_delay_m <= 0) {
        continue;  // Ignore measurements with zero or negative total delay
      }

      m.total_carrier_delay_m
          = dist_bend_m - m.tec_delay_m - second_delay_m / 2.0 - third_delay_m / 3.0;

      // Ephemeris range error (computed later)
      double ephem_range_true
          = (m.posvel_rx_ecef.head<3>() - m.posvel_tx_ecef.head<3>()).norm() - m.clockbias_tx;
      double ephem_range_ephem
          = (m.posvel_rx_ecef.head<3>() - m.posvel_tx_ephem_ecef.head<3>()).norm()
            - m.clockbias_ephem;
      m.ephem_range_error_m = ephem_range_ephem - ephem_range_true;

      // Remove measurements with ephemeris range error is too large (> max_ephem_error_m m)
      if (std::abs(m.ephem_range_error_m) > max_ephem_error_m) {
        continue;
      }

      data.push_back(std::move(m));
      data_idx++;
    }

    // Save the entire data to cache
    save_csv_data_to_cache_(data, cache_file);

    return data;
  }

  void GnssMeasLoader::load(const std::string& csv_path, std::filesystem::path cache_data_path,
                            bool recompute, bool correct_pco, bool correct_clock_bias,
                            double max_ephem_error_m) {
    // Create H5 cache file
    H5Easy::File cache_data_file = GetH5File(cache_data_path, recompute);

    // Save the data to cache
    Container data;

    if (!recompute && cache_data_file.exist("/csv_data")) {
      // Load data from H5 cache
      data = load_h5_data_(cache_data_file);
    } else {
      // For the first time loading from CSV, parse and cache the data
      data = load_csv_data_(csv_path, cache_data_file, recompute, correct_pco, correct_clock_bias,
                            max_ephem_error_m);
    }

    // Convert frames and compute delays
    // Here load from cache if available, other wise, compute frame conversions
    if (!recompute && cache_data_file.exist("/posvel_rx_mci")) {
      std::cout << "  Cache File exists. Loading converted frames from cache..." << std::endl;
      data = load_rv_cache_(data, cache_data_file);
    } else {
      std::cout << "  Cache File does not exist. Computing frame conversions..." << std::endl;
      data = convert_frames_(data, cache_data_file);
    }

    // if data_ is empty, use move assignment (O(1))
    std::cout << "  Assigning loaded GNSS measurements to internal container..." << std::endl;
    if (data_.empty()) {
      data_ = std::move(data);
    } else {
      // Pre-allocate memory to avoid multiple reallocations
      data_.reserve(data_.size() + data.size());

      // Move elements instead of copying them
      data_.insert(data_.end(), std::make_move_iterator(data.begin()),
                   std::make_move_iterator(data.end()));
    }
  }

  GnssMeasLoader::Container GnssMeasLoader::convert_frames_(Container& data,
                                                            H5Easy::File& cache_file) {
    // Converting frames
    std::cout << "  Converting frames for GNSS measurements..." << std::endl;
    int data_count = 0;

    // Load from Cache if available
    auto pbar = Logger::GetProgressBar(data.size(), "Running Frame Conversions", "Main");
    std::atomic<std::size_t> counter{0};

#pragma omp parallel for
    for (std::size_t i = 0; i < data.size(); ++i) {
      auto& m = data[i];

      int tidx = m.tidx;

      // Get Recevier history from pre-loaded matrices
      m.posvel_rx_ecef = posvel_rx_ecef_.row(tidx);
      m.posvel_rx_mci = posvel_rx_mci_.row(tidx);
      m.posvel_rx_pa = posvel_rx_pa_.row(tidx);
      m.posvel_rx_gcrf = posvel_rx_gcrf_.row(tidx);

      // Convert receiver from ECEF to MCI        // Convert to MCI and PA frame later as needed
      m.posvel_tx_mci = ConvertFrame(m.t_tai, m.posvel_tx_ecef, Frame::ECEF, Frame::MOON_CI, false);
      m.posvel_tx_ephem_mci
          = ConvertFrame(m.t_tai, m.posvel_tx_ephem_ecef, Frame::ECEF, Frame::MOON_CI, false);

      // Convert receiver from MCI to PA
      m.posvel_tx_pa
          = ConvertFrame(m.t_tai, m.posvel_tx_mci, Frame::MOON_CI, Frame::MOON_PA, false);
      m.posvel_tx_ephem_pa
          = ConvertFrame(m.t_tai, m.posvel_tx_ephem_mci, Frame::MOON_CI, Frame::MOON_PA, false);

      // Compute GCRF positions
      m.posvel_tx_gcrf = ConvertFrame(m.t_tai, m.posvel_tx_mci, Frame::MOON_CI, Frame::GCRF, false);
      m.posvel_tx_ephem_gcrf
          = ConvertFrame(m.t_tai, m.posvel_tx_ephem_mci, Frame::MOON_CI, Frame::GCRF, false);

      // Compute shapiro delay
      double shapiro_delay = ComputeShapiroDelay(m.t_tai, m.posvel_tx_gcrf.head<3>(),
                                                 m.posvel_rx_gcrf.head<3>(), Frame::GCRF)
                                 .val();
      m.shapiro_delay = shapiro_delay * C;  // convert to meters

      // Compute relativistic delay
      double relativistic_delay
          = ComputeRelativisticDelayLT(m.posvel_rx_mci, Frame::MOON_CI, m.tspan).val();
      m.relativistic_delay = relativistic_delay * C;  // convert to meters

      // Increment atomically
      auto count = counter.fetch_add(1, std::memory_order_relaxed);

      // Let only ONE thread update occasionally (e.g., every 1000 iterations)
      if (count % 10000 == 0) {
#pragma omp critical(progress_update)
        pbar->Update(count);
      }
    }
    pbar->Finish();

    // Save to cache
    save_rv_to_cache_(data, cache_file);

    return data;
  }

  GnssMeasLoader::Container GnssMeasLoader::load_rv_cache_(Container& data,
                                                           H5Easy::File& cache_file) {
    bool recompute = false;

    // Dummy lambda for compatibility
    auto func = []() { return MatX6(); };
    auto func_vec = []() { return VecX(); };

    auto cached_posvel_rx_mci
        = LoadOrRecompute<-1, 6, double>("/posvel_rx_mci", cache_file, recompute, func);
    auto cached_posvel_rx_pa
        = LoadOrRecompute<-1, 6, double>("/posvel_rx_pa", cache_file, recompute, func);
    auto cached_posvel_rx_gcrf
        = LoadOrRecompute<-1, 6, double>("/posvel_rx_gcrf", cache_file, recompute, func);

    auto cached_posvel_tx_mci
        = LoadOrRecompute<-1, 6, double>("/posvel_tx_mci", cache_file, recompute, func);
    auto cached_posvel_tx_pa
        = LoadOrRecompute<-1, 6, double>("/posvel_tx_pa", cache_file, recompute, func);
    auto cached_posvel_tx_gcrf
        = LoadOrRecompute<-1, 6, double>("/posvel_tx_gcrf", cache_file, recompute, func);
    // Ephemeris-based
    auto cached_posvel_tx_ephem_mci
        = LoadOrRecompute<-1, 6, double>("/posvel_tx_ephem_mci", cache_file, recompute, func);
    auto cached_posvel_tx_ephem_gcrf
        = LoadOrRecompute<-1, 6, double>("/posvel_tx_ephem_gcrf", cache_file, recompute, func);
    auto cached_posvel_tx_ephem_pa
        = LoadOrRecompute<-1, 6, double>("/posvel_tx_ephem_pa", cache_file, recompute, func);
    // Delays
    auto cached_shapiro_delays
        = LoadOrRecompute<-1, 1, double>("/shapiro_delays", cache_file, recompute, func_vec);
    auto cached_relativistic_delays
        = LoadOrRecompute<-1, 1, double>("/relativistic_delays", cache_file, recompute, func_vec);

    // throw error if its recomputing here (should not happen)
    if (cached_posvel_tx_mci.rows() != data.size()) {
      throw std::runtime_error("Cached converted frame size mismatch.");
    }

    for (size_t i = 0; i < data.size(); ++i) {
      int tidx = data[i].tidx;
      // For receiver, use tidx to index
      data[i].posvel_rx_mci = cached_posvel_rx_mci.row(tidx).transpose();
      data[i].posvel_rx_pa = cached_posvel_rx_pa.row(tidx).transpose();
      data[i].posvel_rx_gcrf = cached_posvel_rx_gcrf.row(tidx).transpose();
      // For transmitter, use data index
      data[i].posvel_tx_mci = cached_posvel_tx_mci.row(static_cast<int>(i)).transpose();
      data[i].posvel_tx_pa = cached_posvel_tx_pa.row(static_cast<int>(i)).transpose();
      data[i].posvel_tx_gcrf = cached_posvel_tx_gcrf.row(static_cast<int>(i)).transpose();
      // Ephemeris-based
      data[i].posvel_tx_ephem_mci = cached_posvel_tx_ephem_mci.row(static_cast<int>(i)).transpose();
      data[i].posvel_tx_ephem_gcrf
          = cached_posvel_tx_ephem_gcrf.row(static_cast<int>(i)).transpose();
      data[i].posvel_tx_ephem_pa = cached_posvel_tx_ephem_pa.row(static_cast<int>(i)).transpose();
      data[i].shapiro_delay = cached_shapiro_delays(i);
      data[i].relativistic_delay = cached_relativistic_delays(i);
    }
    std::cout << "  Loaded converted frames from cache." << std::endl;

    return data;
  }

  void GnssMeasLoader::save_rv_to_cache_(Container& data, H5Easy::File& cache_file) {
    // Save histories to cache
    int data_idx = 0;

    construct_posvel_histories_(data);

    // Receiver history
    Dump(cache_file, "/posvel_rx_ecef", posvel_rx_ecef_, H5Easy::DumpMode::Overwrite);
    Dump(cache_file, "/posvel_rx_mci", posvel_rx_mci_, H5Easy::DumpMode::Overwrite);
    Dump(cache_file, "/posvel_rx_pa", posvel_rx_pa_, H5Easy::DumpMode::Overwrite);
    Dump(cache_file, "/posvel_rx_gcrf", posvel_rx_gcrf_, H5Easy::DumpMode::Overwrite);
    // True transmitter history
    Dump(cache_file, "/posvel_tx_ecef", posvel_tx_ecef_, H5Easy::DumpMode::Overwrite);
    Dump(cache_file, "/posvel_tx_mci", posvel_tx_mci_, H5Easy::DumpMode::Overwrite);
    Dump(cache_file, "/posvel_tx_pa", posvel_tx_pa_, H5Easy::DumpMode::Overwrite);
    Dump(cache_file, "/posvel_tx_gcrf", posvel_tx_gcrf_, H5Easy::DumpMode::Overwrite);
    // Ephemeris-based
    Dump(cache_file, "/posvel_tx_ephem_ecef", posvel_tx_ephem_ecef_, H5Easy::DumpMode::Overwrite);
    Dump(cache_file, "/posvel_tx_ephem_mci", posvel_tx_ephem_mci_, H5Easy::DumpMode ::Overwrite);
    Dump(cache_file, "/posvel_tx_ephem_gcrf", posvel_tx_ephem_gcrf_, H5Easy::DumpMode::Overwrite);
    Dump(cache_file, "/posvel_tx_ephem_pa", posvel_tx_ephem_pa_, H5Easy::DumpMode::Overwrite);
    // Shapiro and relativistic delays
    Dump(cache_file, "/shapiro_delays", shapiro_delays_, H5Easy::DumpMode::Overwrite);
    Dump(cache_file, "/relativistic_delays", relativistic_delays_, H5Easy::DumpMode::Overwrite);
  }

  void GnssMeasLoader::save_csv_data_to_cache_(Container& data, H5Easy::File& cache_file) {
    int num_cols = 37;

    MatX datamat(data.size(), num_cols);

    std::cout << "   Saving " << data.size() << " GNSS measurements to cache." << std::endl;

    for (size_t i = 0; i < data.size(); ++i) {
      // Fill datamat with CSV data from each measurement
      // Assuming each measurement has a method to_csv_row() returning Eigen::RowVectorXd
      VecX data_row(num_cols);

      data_row(0) = data[i].tidx;
      data_row(1) = data[i].t_tai;
      data_row(2) = data[i].tspan;
      data_row(3) = data[i].epoch_t;
      data_row(4) = gnss_const_to_double(data[i].gnss_const);
      data_row(5) = signal_to_double(data[i].signal);
      data_row(6) = data[i].sat_id;
      data_row(7) = data[i].prn;
      data_row(8) = data[i].min_alt;
      data_row(9) = data[i].sigma_range;
      data_row(10) = data[i].sigma_rangerate;
      data_row(11) = data[i].sigma_carrier;
      data_row(12) = data[i].cn0;
      data_row.segment<3>(13) = data[i].posvel_tx_ecef.segment<3>(0);
      data_row.segment<3>(16) = data[i].posvel_rx_ecef.segment<3>(0);
      data_row.segment<3>(19) = data[i].posvel_tx_ecef.segment<3>(3);
      data_row.segment<3>(22) = data[i].posvel_rx_ecef.segment<3>(3);
      data_row(25) = data[i].clockbias_tx;
      data_row.segment<3>(26) = data[i].posvel_tx_ephem_ecef.segment<3>(0);
      data_row.segment<3>(29) = data[i].posvel_tx_ephem_ecef.segment<3>(3);
      data_row(32) = data[i].clockbias_ephem;
      data_row(33) = data[i].tec_delay_m;
      data_row(34) = data[i].total_delay_m;
      data_row(35) = data[i].total_carrier_delay_m;
      data_row(36) = data[i].ephem_range_error_m;
      datamat.row(static_cast<int>(i)) = data_row;
    }
    // Save CSV data to cache
    Dump(cache_file, "/csv_data", datamat, H5Easy::DumpMode::Overwrite);
  }

  GnssMeasLoader::Container GnssMeasLoader::load_h5_data_(H5Easy::File& cache_file) {
    bool recompute = false;
    MatX datamat;

    auto func = []() { return MatX(); };

    auto cached_csv_data
        = LoadOrRecompute<-1, 37, double>("/csv_data", cache_file, recompute, func);

    int data_size = static_cast<int>(cached_csv_data.rows());

    Container data(data_size);

    std::cout << "   Loading " << data_size << " GNSS measurements from cache." << std::endl;

    for (int i = 0; i < data_size; ++i) {
      GnssMeas m;
      m.tidx = static_cast<int>(cached_csv_data(i, 0));
      m.t_tai = cached_csv_data(i, 1);
      m.tspan = cached_csv_data(i, 2);
      m.epoch_t = cached_csv_data(i, 3);
      m.gnss_const = gnss_const_to_string(cached_csv_data(i, 4));
      m.signal = signal_to_string(cached_csv_data(i, 5));
      m.sat_id = static_cast<int>(cached_csv_data(i, 6));
      m.prn = static_cast<int>(cached_csv_data(i, 7));
      m.min_alt = cached_csv_data(i, 8);
      m.sigma_range = cached_csv_data(i, 9);
      m.sigma_rangerate = cached_csv_data(i, 10);
      m.sigma_carrier = cached_csv_data(i, 11);
      m.cn0 = cached_csv_data(i, 12);
      m.posvel_tx_ecef.segment<3>(0) = cached_csv_data.row(i).segment<3>(13);
      m.posvel_rx_ecef.segment<3>(0) = cached_csv_data.row(i).segment<3>(16);
      m.posvel_tx_ecef.segment<3>(3) = cached_csv_data.row(i).segment<3>(19);
      m.posvel_rx_ecef.segment<3>(3) = cached_csv_data.row(i).segment<3>(22);
      m.clockbias_tx = cached_csv_data(i, 25);
      m.posvel_tx_ephem_ecef.segment<3>(0) = cached_csv_data.row(i).segment<3>(26);
      m.posvel_tx_ephem_ecef.segment<3>(3) = cached_csv_data.row(i).segment<3>(29);
      m.clockbias_ephem = cached_csv_data(i, 32);
      m.tec_delay_m = cached_csv_data(i, 33);
      m.total_delay_m = cached_csv_data(i, 34);
      m.total_carrier_delay_m = cached_csv_data(i, 35);
      m.ephem_range_error_m = cached_csv_data(i, 36);
      data[i] = m;
    }

    return data;
  }

  void GnssMeasLoader::clear_history_vectors() {
    tidx_to_dataidx_.clear();
    tidx_to_const_signal_prn_.clear();

    // matrices
    posvel_tx_ecef_.resize(0, 0);
    posvel_tx_mci_.resize(0, 0);
    posvel_tx_pa_.resize(0, 0);
    posvel_tx_gcrf_.resize(0, 0);

    posvel_tx_ephem_ecef_.resize(0, 0);
    posvel_tx_ephem_mci_.resize(0, 0);
    posvel_tx_ephem_gcrf_.resize(0, 0);
    posvel_tx_ephem_pa_.resize(0, 0);

    shapiro_delays_.resize(0);
    relativistic_delays_.resize(0);
  }

  // After reading all lines, convert vectors to Eigen types --------------------------------
  void GnssMeasLoader::construct_posvel_histories_(Container& data) {
    clear_history_vectors();  // Do not accumulate

    // resize posvel matrices
    posvel_tx_ecef_.resize(data.size(), 6);
    posvel_tx_mci_.resize(data.size(), 6);
    posvel_tx_pa_.resize(data.size(), 6);
    posvel_tx_gcrf_.resize(data.size(), 6);
    posvel_tx_ephem_ecef_.resize(data.size(), 6);
    posvel_tx_ephem_mci_.resize(data.size(), 6);
    posvel_tx_ephem_gcrf_.resize(data.size(), 6);
    posvel_tx_ephem_pa_.resize(data.size(), 6);
    shapiro_delays_.resize(data.size());
    relativistic_delays_.resize(data.size());

    int i = 0;
    for (const auto& m : data) {
      tidx_to_dataidx_[m.tidx].push_back(i);
      tidx_to_const_signal_prn_[m.tidx].emplace_back(
          ConstFreqPrnTuple{m.gnss_const, m.signal, m.prn});

      // Always store Tx history
      posvel_tx_ecef_.row(static_cast<int>(i)) = m.posvel_tx_ecef.transpose();
      posvel_tx_mci_.row(static_cast<int>(i)) = m.posvel_tx_mci.transpose();
      posvel_tx_pa_.row(static_cast<int>(i)) = m.posvel_tx_pa.transpose();
      posvel_tx_gcrf_.row(static_cast<int>(i)) = m.posvel_tx_gcrf.transpose();

      posvel_tx_ephem_ecef_.row(static_cast<int>(i)) = m.posvel_tx_ephem_ecef.transpose();
      posvel_tx_ephem_mci_.row(static_cast<int>(i)) = m.posvel_tx_ephem_mci.transpose();
      posvel_tx_ephem_gcrf_.row(static_cast<int>(i)) = m.posvel_tx_ephem_gcrf.transpose();
      posvel_tx_ephem_pa_.row(static_cast<int>(i)) = m.posvel_tx_ephem_pa.transpose();
      shapiro_delays_(i) = m.shapiro_delay;
      relativistic_delays_(i) = m.relativistic_delay;
      i++;
    }
  }

  void GnssMeasLoader::print_delay_statistics() {
    std::cout << " " << std::endl;
    std::cout << "[GNSS Measurement Data Statistics]" << std::endl;
    std::cout << " GNSS Measurements Loaded: " << data_.size() << " total measurements."
              << std::endl;

    int data_size = static_cast<int>(data_.size());
    std::vector<double> shapiro_delays(data_size);
    std::vector<double> relativistic_delays(data_size);
    std::vector<double> iono_delays(data_size);
    std::vector<double> ephem_range_errors(data_size);
    std::vector<double> sigma_ranges(data_size);
    std::vector<double> sigma_carriers(data_size);

    for (int i = 0; i < data_size; ++i) {
      const auto& m = data_[i];
      shapiro_delays[i] = m.shapiro_delay;
      relativistic_delays[i] = m.relativistic_delay;
      iono_delays[i] = m.total_delay_m;
      ephem_range_errors[i] = m.ephem_range_error_m;
      sigma_ranges[i] = m.sigma_range;
      sigma_carriers[i] = m.sigma_carrier;
    }

    // Compute Mean, std dev, median, 95%, 99%
    auto compute_stats = [](const std::vector<double>& vec) {
      // Mean
      double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
      double mean = sum / static_cast<double>(vec.size());
      // Std Dev
      double sq_sum = std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
      double std_dev = std::sqrt(sq_sum / static_cast<double>(vec.size()) - mean * mean);
      // For median and percentiles, we need to sort the abs values
      std::vector<double> abs_vec(vec.size());
      std::transform(vec.begin(), vec.end(), abs_vec.begin(), [](double x) { return std::abs(x); });
      // Median
      std::vector<double> sorted_vec = abs_vec;  // copy for sorting
      std::sort(sorted_vec.begin(), sorted_vec.end());
      double median;
      size_t n = sorted_vec.size();
      if (n % 2 == 0) {
        median = 0.5 * (sorted_vec[n / 2 - 1] + sorted_vec[n / 2]);
      } else {
        median = sorted_vec[n / 2];
      }

      // 95th percentile
      double p95 = sorted_vec[static_cast<size_t>(0.95 * n)];
      // 99th percentile
      double p99 = sorted_vec[static_cast<size_t>(0.99 * n)];

      return std::make_tuple(mean, std_dev, median, p95, p99);
    };

    auto [shapiro_mean, shapiro_std, shapiro_median, shapiro_p95, shapiro_p99]
        = compute_stats(shapiro_delays);
    auto [relativistic_mean, relativistic_std, relativistic_median, relativistic_p95,
          relativistic_p99]
        = compute_stats(relativistic_delays);
    auto [iono_mean, iono_std, iono_median, iono_p95, iono_p99] = compute_stats(iono_delays);
    auto [ephem_mean, ephem_std, ephem_median, ephem_p95, ephem_p99]
        = compute_stats(ephem_range_errors);
    auto [sigma_range_mean, sigma_range_std, sigma_range_median, sigma_range_p95, sigma_range_p99]
        = compute_stats(sigma_ranges);
    auto [sigma_carrier_mean, sigma_carrier_std, sigma_carrier_median, sigma_carrier_p95,
          sigma_carrier_p99]
        = compute_stats(sigma_carriers);

    // Print statistics as tables
    std::cout << " ---------------------------------------------------------------------"
              << std::endl;
    std::cout << " Delay Type          | Mean (m)      | Std Dev (m)  | Median (m)  | 95th % (m)  "
                 "| 99th % (m)"
              << std::endl;
    std::cout << "---------------------------------------------------------------------"
              << std::endl;
    std::cout << " Shapiro Delay       | " << shapiro_mean << " | " << shapiro_std << " | "
              << shapiro_median << " | " << shapiro_p95 << " | " << shapiro_p99 << std::endl;
    std::cout << " Relativistic Delay  | " << relativistic_mean << " | " << relativistic_std
              << " | " << relativistic_median << " | " << relativistic_p95 << " | "
              << relativistic_p99 << std::endl;
    std::cout << " Ionospheric Delay   | " << iono_mean << " | " << iono_std << " | " << iono_median
              << " | " << iono_p95 << " | " << iono_p99 << std::endl;
    std::cout << " Ephemeris Range Err | " << ephem_mean << " | " << ephem_std << " | "
              << ephem_median << " | " << ephem_p95 << " | " << ephem_p99 << std::endl;
    std::cout << " Sigma Range         | " << sigma_range_mean << " | " << sigma_range_std << " | "
              << sigma_range_median << " | " << sigma_range_p95 << " | " << sigma_range_p99
              << std::endl;
    std::cout << " Sigma Carrier       | " << sigma_carrier_mean << " | " << sigma_carrier_std
              << " | " << sigma_carrier_median << " | " << sigma_carrier_p95 << " | "
              << sigma_carrier_p99 << std::endl;
    std::cout << "---------------------------------------------------------------------"
              << std::endl;
    std::cout << " " << std::endl;
  }

  GnssMeasLoader::Container GnssMeasLoader::filter_by_prn(const Container& all, int prn) {
    Container out;
    out.reserve(all.size() / 4);  // rough guess; not critical

    for (const auto& m : all) {
      if (m.prn == prn) {
        out.push_back(m);
      }
    }
    out.shrink_to_fit();
    return out;
  }

  GnssMeasLoader::TidxGroups GnssMeasLoader::group_by_tidx(const Container& all) {
    TidxGroups groups;

    for (const auto& m : all) {
      groups[m.tidx].push_back(m);
    }

    return groups;
  }

}  // namespace filtering_sim
