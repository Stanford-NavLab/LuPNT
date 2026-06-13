#include "lupnt/data/kernels.h"

#include <filesystem>
#include <mutex>
#include <string>

#include "lupnt/conversions/frame_converter.h"
#include "lupnt/conversions/time_conversions.h"
#include "lupnt/data/kernels.h"

namespace lupnt {

  namespace EphemID {
    int MERCURY_BARYCENTER = 0;
    int VENUS_BARYCENTER = 1;
    int EARTH_MOON_BARYCENTER = 2;
    int EMB = EARTH_MOON_BARYCENTER;
    int MARS_BARYCENTER = 3;
    int JUPITER_BARYCENTER = 4;
    int SATURN_BARYCENTER = 5;
    int URANUS_BARYCENTER = 6;
    int NEPTUNE_BARYCENTER = 7;
    int PLUTO_BARYCENTER = 8;
    int MOON = 9;
    int SUN = 10;
    int EARTH_NUTATIONS = 11;
    int MOON_MANTLE_LIBRATIONS = 12;
    int MOON_MANTLE_ANGULAR_VELOCITY = 13;
    int TT_TDB = 14;
  };  // namespace EphemID

  std::map<BodyId, int> naif2ephemId = {
      {BodyId::MERCURY, EphemID::MERCURY_BARYCENTER},
      {BodyId::MERCURY_BARYCENTER, EphemID::MERCURY_BARYCENTER},
      {BodyId::VENUS, EphemID::VENUS_BARYCENTER},
      {BodyId::VENUS_BARYCENTER, EphemID::VENUS_BARYCENTER},
      {BodyId::MARS, EphemID::MARS_BARYCENTER},
      {BodyId::MARS_BARYCENTER, EphemID::MARS_BARYCENTER},
      {BodyId::EARTH_MOON_BARYCENTER, EphemID::EARTH_MOON_BARYCENTER},
      {BodyId::JUPITER, EphemID::JUPITER_BARYCENTER},
      {BodyId::JUPITER_BARYCENTER, EphemID::JUPITER_BARYCENTER},
      {BodyId::SATURN, EphemID::SATURN_BARYCENTER},
      {BodyId::SATURN_BARYCENTER, EphemID::SATURN_BARYCENTER},
      {BodyId::URANUS, EphemID::URANUS_BARYCENTER},
      {BodyId::URANUS_BARYCENTER, EphemID::URANUS_BARYCENTER},
      {BodyId::NEPTUNE, EphemID::NEPTUNE_BARYCENTER},
      {BodyId::NEPTUNE_BARYCENTER, EphemID::NEPTUNE_BARYCENTER},
      {BodyId::PLUTO_BARYCENTER, EphemID::PLUTO_BARYCENTER},
      {BodyId::MOON, EphemID::MOON},
      {BodyId::SUN, EphemID::SUN},
  };

  namespace {
    Vec6 StateFromSI(const Vec6& rv, const UnitSystem& units) {
      Vec6 out;
      out.head(3) = rv.head(3) / units.length;
      out.tail(3) = rv.tail(3) * (units.time / units.length);
      return out;
    }

    Vec3 PositionFromSI(const Vec3& r, const UnitSystem& units) { return r / units.length; }
  }  // namespace

  struct EphemerisHeaderData {
    int KSIZE;
    int NCOEFF;
    std::string info;
    std::string start_info;
    std::string final_info;
    double jd_tdb_start;
    double jd_tdb_end;
    double step;
    std::vector<int> coeff_offset;
    std::vector<int> n_coeffs;
    std::vector<int> n_subintervals;
    std::map<std::string, double> constants;
    std::vector<int> n_properties = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 1};
    double emrat;
  };

  struct EphemerisBlock {
    double jd_tdb_start;
    double jd_tdb_end;
    double coeff[MAXCOEFF];
  };

  struct EphemerisData {
    EphemerisHeaderData header;
    double jd_tdb_start;
    double jd_tdb_end;
    std::vector<EphemerisBlock> blocks;
  };

  UniquePtr<EphemerisData> ephemeris_data;
  std::mutex kernels_mutex;

  double ParseDouble(const std::string& str) {
    std::string s = str;
    for (char& c : s)
      if (c == 'D') c = 'e';
    return std::stod(s);
  }

  void ParseGroup1010(std::ifstream& infile, EphemerisHeaderData& data) {
    std::string line, empty_line;
    std::getline(infile, empty_line);
    std::getline(infile, data.info);
    std::getline(infile, data.start_info);
    std::getline(infile, data.final_info);
    std::getline(infile, empty_line);
  }

  void ParseGroup1030(std::ifstream& infile, EphemerisHeaderData& data) {
    std::string line, empty_line;
    std::getline(infile, empty_line);
    std::getline(infile, line);
    std::istringstream iss(line);
    iss >> data.jd_tdb_start >> data.jd_tdb_end >> data.step;
    std::getline(infile, empty_line);
  }

  std::vector<std::string> ParseGroup1040(std::ifstream& infile, EphemerisHeaderData& data) {
    (void)data;
    std::string line, empty_line;
    std::getline(infile, empty_line);
    int num_constants;
    infile >> num_constants;
    std::string constant;
    std::vector<std::string> constant_names;
    for (int i = 0; i < num_constants; ++i) {
      infile >> constant;
      constant_names.push_back(constant);
    }
    std::getline(infile, empty_line);
    std::getline(infile, empty_line);
    return constant_names;
  }

  std::vector<double> ParseGroup1041(std::ifstream& infile, EphemerisHeaderData& data) {
    (void)data;
    std::string line, empty_line;
    std::getline(infile, empty_line);
    int num_values;
    infile >> num_values;
    std::string value_str;
    std::vector<double> constant_values;
    for (int i = 0; i < num_values; ++i) {
      infile >> value_str;
      constant_values.push_back(ParseDouble(value_str));
    }
    std::getline(infile, empty_line);
    std::getline(infile, empty_line);
    return constant_values;
  }

  void ParseGroup1050(std::ifstream& infile, EphemerisHeaderData& data) {
    std::string line, empty_line;
    std::getline(infile, empty_line);
    std::vector<int> start_locations, n_coeffs, n_subintervals;
    // Read the next three lines and store the values in respective vectors
    for (int i = 0; i < 3; ++i) {
      std::getline(infile, line);
      std::istringstream iss(line);
      int value;
      while (iss >> value) {
        if (i == 0) {
          start_locations.push_back(value);
        } else if (i == 1) {
          n_coeffs.push_back(value);
        } else if (i == 2) {
          n_subintervals.push_back(value);
        }
      }
    }
    std::getline(infile, empty_line);
    data.coeff_offset = start_locations;
    data.n_coeffs = n_coeffs;
    data.n_subintervals = n_subintervals;
  }

  void ReadEphemerisHeaderFile(const std::filesystem::path& filepath, EphemerisHeaderData& data) {
    std::ifstream infile(filepath);
    if (!infile.is_open()) {
      std::string msg = "Unable to open file: " + filepath.string();
      throw std::runtime_error(msg);
    }
    std::string line, empty_line;
    std::vector<std::string> constant_names;
    std::vector<double> constant_values;
    while (std::getline(infile, line)) {
      if (line.find("KSIZE") != std::string::npos) {
        std::istringstream iss(line);
        std::string temp;
        iss >> temp >> data.KSIZE >> temp >> data.NCOEFF;
      } else if (line.find("GROUP   1010") != std::string::npos) {
        ParseGroup1010(infile, data);
      } else if (line.find("GROUP   1030") != std::string::npos) {
        ParseGroup1030(infile, data);
      } else if (line.find("GROUP   1040") != std::string::npos) {
        constant_names = ParseGroup1040(infile, data);
      } else if (line.find("GROUP   1041") != std::string::npos) {
        constant_values = ParseGroup1041(infile, data);
      } else if (line.find("GROUP   1050") != std::string::npos) {
        ParseGroup1050(infile, data);
      }
    }
    infile.close();
    for (size_t i = 0; i < constant_names.size(); ++i)
      data.constants[constant_names[i]] = constant_values[i];
    data.emrat = data.constants["EMRAT"];
  }

  // void ReadEphemerisCoefficientsFile(const std::filesystem::path& filepath,
  //                                    EphemerisHeaderData& data) {
  //   ephemeris_data = MakeUnique<EphemerisData>();
  //   ephemeris_data->header = data;
  //   ephemeris_data->blocks.clear();

  //   std::ifstream infile(filepath);
  //   if (!infile.is_open()) {
  //     std::string msg = "Unable to open file: " + filepath.string();
  //     throw std::runtime_error(msg);
  //   }
  //   std::cout << "Loading ephemeris coefficients from " + filepath.string() + " ..." <<
  //   std::endl; std::string value_str, line; int block = 1; int block_in, tmp; while
  //   (std::getline(infile, line)) {
  //     std::istringstream iss(line);
  //     iss >> block_in >> tmp;
  //     std::cout << "\r  Reading block " << block_in << std::flush;
  //     std::cout << "Block expected: " << block << ", Block read: " << block_in << std::endl;
  //     if (block_in != block) throw std::runtime_error("Block number mismatch");
  //     EphemerisBlock block_data;
  //     infile >> value_str;
  //     block_data.jd_tdb_start = ParseDouble(value_str);
  //     infile >> value_str;
  //     block_data.jd_tdb_end = ParseDouble(value_str);
  //     for (int i = 0; i < data.NCOEFF - 2; ++i) {
  //       infile >> value_str;
  //       block_data.coeff[i] = ParseDouble(value_str);
  //     }
  //     ephemeris_data->blocks.push_back(block_data);
  //     std::getline(infile, line);  // Finish the line
  //     block++;
  //   }
  //   infile.close();
  //   ephemeris_data->jd_tdb_start = ephemeris_data->blocks[0].jd_tdb_start;
  //   ephemeris_data->jd_tdb_end = ephemeris_data->blocks.back().jd_tdb_end;
  // }

  void ReadEphemerisCoefficientsFile(const std::filesystem::path& filepath,
                                     EphemerisHeaderData& data) {
    ephemeris_data = MakeUnique<EphemerisData>();
    ephemeris_data->header = data;
    ephemeris_data->blocks.clear();

    std::ifstream infile(filepath);
    if (!infile.is_open()) {
      std::string msg = "Unable to open file: " + filepath.string();
      throw std::runtime_error(msg);
    }
    // std::cout << "Loading ephemeris coefficients from "
    //           << filepath.string() << " ..." << std::endl;

    std::string value_str;
    int block = 1;
    int block_in, tmp;

    // Read until we fail to read a new block header
    while (infile >> block_in >> tmp) {
      // std::cout << "  Reading block " << block_in
      //           << " (expected " << block << ")\n";

      // Basic sanity check on block numbering
      if (block_in != block) {
        throw std::runtime_error("Block number mismatch");
      }

      EphemerisBlock block_data;

      // Start / end JD
      infile >> value_str;
      block_data.jd_tdb_start = ParseDouble(value_str);

      infile >> value_str;
      block_data.jd_tdb_end = ParseDouble(value_str);

      // Coefficients (NCOEFF includes the 2 JD entries)
      for (int i = 0; i < data.NCOEFF - 2; ++i) {
        infile >> value_str;
        block_data.coeff[i] = ParseDouble(value_str);
      }

      ephemeris_data->blocks.push_back(block_data);
      ++block;
    }

    infile.close();

    if (ephemeris_data->blocks.empty()) {
      throw std::runtime_error("No ephemeris blocks read from file: " + filepath.string());
    }

    ephemeris_data->jd_tdb_start = ephemeris_data->blocks.front().jd_tdb_start;
    ephemeris_data->jd_tdb_end = ephemeris_data->blocks.back().jd_tdb_end;
  }

  Vec2 ComputePolynomialPosVel(Real x, const double* scale, const double* coeff, int offset,
                               int num) {
    Real x2, w0 = 0., w1 = 0., dw0 = 0., dw1 = 0., tmp;

    x = (x - scale[0]) / scale[1];
    x2 = x * 2.;
    while (--num) {
      tmp = dw1;
      dw1 = dw0;
      dw0 = w0 * 2. + dw0 * x2 - tmp;
      tmp = w1;
      w1 = w0;
      w0 = coeff[offset + num] + (x2 * w0 - tmp);
    }
    Real f = coeff[offset] + (x * w0 - w1);
    Real df = (w0 + x * dw0 - dw1) / scale[1];
    return Vec2(f, df) * M_KM;
  }

  Real ComputePolynomialPos(Real x, const double* scale, const double* coeff, int offset, int num) {
    Real x2, w0 = 0., w1 = 0., tmp;

    x = (x - scale[0]) / scale[1];
    x2 = x * 2.;
    while (--num) {
      tmp = w1;
      w1 = w0;
      w0 = coeff[offset + num] + (x2 * w0 - tmp);
    }
    Real f = coeff[offset] + (x * w0 - w1);
    return f * M_KM;
  }

  void LoadEphemerisData() {
    std::lock_guard<std::mutex> lock(kernels_mutex);
    if (ephemeris_data) return;  // Data already loaded

    EphemerisHeaderData data;
    ReadEphemerisHeaderFile(GetAsciiKernelDir() / "de440" / "header.440", data);
    ReadEphemerisCoefficientsFile(GetAsciiKernelDir() / "de440" / "ascp01950.440", data);
  }

  struct CacheKernelData {
    Real t_tdb = NAN;
    Vec6 rv;
    bool compute_vel;
  };
  static std::array<CacheKernelData, 15> cache_kernel_data;
  static std::array<std::mutex, 15> cache_mutex;

  Vec6 GetKernelData(Real t_tdb, int id, bool compute_vel = true) {
    if (!ephemeris_data) LoadEphemerisData();

    // Check cache
    // std::unique_lock<std::mutex> lock(cache_mutex[id]);
    // CacheKernelData cache = cache_kernel_data[id];
    // if (abs(cache.t_tdb - t_tdb) < EPS && cache.compute_vel == compute_vel) {
    //   lock.unlock();
    //   return cache.rv;
    // }
    // lock.unlock();

    Real jd_tdb = TimeToJd(t_tdb);

    // Block
    if (jd_tdb < ephemeris_data->jd_tdb_start || jd_tdb > ephemeris_data->jd_tdb_end) {
      std::string msg = "Kernels not loaded for the requested time\n";
      msg += "- Start: " + std::to_string(ephemeris_data->jd_tdb_start) + "\n";
      msg += "- End: " + std::to_string(ephemeris_data->jd_tdb_end) + "\n";
      msg += "- Requested: " + std::to_string(jd_tdb.val());
      throw std::runtime_error(msg);
    }
    double Dt = ephemeris_data->header.step;
    int i = int((jd_tdb - ephemeris_data->jd_tdb_start) / Dt);
    if (!(i >= 0 && i < (int)ephemeris_data->blocks.size())) {
      std::string msg = "Block index out of range\n";
      msg += "- Index: " + std::to_string(i) + "\n";
      msg += "- Number of blocks: " + std::to_string(ephemeris_data->blocks.size());
      throw std::runtime_error(msg);
    }
    EphemerisBlock& block = ephemeris_data->blocks[i];
    EphemerisHeaderData& header = ephemeris_data->header;

    // Subinterval
    double Dt_subint = Dt / header.n_subintervals[id];
    int j = int((jd_tdb - block.jd_tdb_start) / Dt_subint);
    int n_coeff = header.n_coeffs[id] * header.n_properties[id];
    int offset = header.coeff_offset[id] + j * n_coeff - 3.;
    double jd_tdb_subint = block.jd_tdb_start + j * Dt_subint;
    double scale[2] = {jd_tdb_subint + Dt_subint / 2., Dt_subint / 2.};  // center and half width
    scale[0] = JdToTime(scale[0]).val();
    scale[1] *= SECS_DAY;

    Vec6 rv;
    for (int i = 0; i < 3; i++) {
      if (compute_vel) {
        Vec2 rv_tmp
            = ComputePolynomialPosVel(t_tdb, scale, block.coeff, offset, header.n_coeffs[id]);
        rv[i] = rv_tmp(0);
        rv[i + 3] = rv_tmp(1);
      } else {
        rv[i] = ComputePolynomialPos(t_tdb, scale, block.coeff, offset, header.n_coeffs[id]);
        rv[i + 3] = 0.;
      }
      offset += header.n_coeffs[id];
    }

    // cache.t_tdb = t_tdb;
    // cache.rv = rv;
    // cache.compute_vel = compute_vel;

    // lock.lock();
    // if (t_tdb > cache_kernel_data[id].t_tdb + EPS) {
    //   cache_kernel_data[id] = cache;
    // }
    // lock.unlock();

    return rv;
  }  // namespace lupnt

  Vec6 GetLunarMantleData(Real t_tdb, bool compute_vel) {
    return GetKernelData(t_tdb, EphemID::MOON_MANTLE_LIBRATIONS, compute_vel) / M_KM;
  }

  MatX6 GetLunarMantleData(VecX t_tdb, bool compute_vel) {
    MatX6 rv(t_tdb.size(), 6);
    for (int i = 0; i < t_tdb.size(); i++) {
      rv.row(i) = GetLunarMantleData(t_tdb(i), compute_vel) / M_KM;
    }
    return rv;
  }

  Vec6 GetBodyPosVel(Real t_tdb, BodyId center, BodyId target, Frame frame) {
    if (!ephemeris_data) LoadEphemerisData();

    if (center == target) return Vec6::Zero();
    Vec6 rv_center = Vec6::Zero();
    Vec6 rv_target = Vec6::Zero();
    double emr = ephemeris_data->header.emrat;

    // Earth-Moon system
    if (center == BodyId::EARTH && target == BodyId::MOON) {
      rv_target = GetKernelData(t_tdb, EphemID::MOON);
    } else if (center == BodyId::MOON && target == BodyId::EARTH) {
      rv_center = GetKernelData(t_tdb, EphemID::MOON);

    } else if (center == BodyId::EMB && target == BodyId::MOON) {
      rv_target = GetKernelData(t_tdb, EphemID::MOON);
      rv_center = rv_target / (1. + emr);
    } else if (center == BodyId::MOON && target == BodyId::EMB) {
      rv_center = GetKernelData(t_tdb, EphemID::MOON);
      rv_target = rv_center / (1. + emr);

    } else if (center == BodyId::EMB && target == BodyId::EARTH) {
      rv_center = GetKernelData(t_tdb, EphemID::MOON) / (1. + emr);
    } else if (center == BodyId::EARTH && target == BodyId::EMB) {
      rv_target = GetKernelData(t_tdb, EphemID::MOON) / (1. + emr);
    } else {
      if (center == BodyId::EARTH || center == BodyId::MOON || target == BodyId::EARTH
          || target == BodyId::MOON) {
        Vec6 rv_emb = GetKernelData(t_tdb, EphemID::EMB);
        Vec6 rv_moon = GetKernelData(t_tdb, EphemID::MOON);
        Vec6 rv_earth = rv_emb - rv_moon / (1. + emr);

        if (center == BodyId::EARTH)
          rv_center = rv_earth;
        else if (center == BodyId::MOON)
          rv_center = rv_earth + rv_moon;

        if (target == BodyId::EARTH)
          rv_target = rv_earth;
        else if (target == BodyId::MOON)
          rv_target = rv_earth + rv_moon;
      }

      if (center != BodyId::EARTH && center != BodyId::MOON && center != BodyId::SSB)
        rv_center = GetKernelData(t_tdb, naif2ephemId.at(center));

      if (target != BodyId::EARTH && target != BodyId::MOON && target != BodyId::SSB)
        rv_target = GetKernelData(t_tdb, naif2ephemId.at(target));
    }

    if (frame == Frame::GCRF) return rv_target - rv_center;
    rv_target = ConvertFrame(t_tdb, rv_target, Frame::GCRF, frame);
    rv_center = ConvertFrame(t_tdb, rv_center, Frame::GCRF, frame);
    return rv_target - rv_center;
  }

  Vec3 GetBodyPos(Real t_tdb, BodyId center, BodyId target, Frame frame) {
    return GetBodyPosVel(t_tdb, center, target, frame).head(3);
  }

  Vec6 GetBodyPosVel(Real t_tdb, BodyId target, Frame frame) {
    LUPNT_CHECK(frame != Frame::UNDEFINED, "Frame is undefined", "GetBodyPosVel");
    return GetBodyPosVel(t_tdb, frame_centers.at(frame), target, frame);
  }

  Vec6 GetBodyPosVel(Real t_tdb, BodyId center, BodyId target, Frame frame, const UnitSystem& units,
                     CoordinateScale scale) {
    Vec6 rv = GetBodyPosVel(t_tdb, center, target, frame);
    rv = ScaleStateForCoordinateScale(rv, CoordinateScale::TDB, scale);
    return StateFromSI(rv, units);
  }

  Vec6 GetBodyPosVel(Real t_tdb, BodyId target, Frame frame, const UnitSystem& units,
                     CoordinateScale scale) {
    LUPNT_CHECK(frame != Frame::UNDEFINED, "Frame is undefined", "GetBodyPosVel");
    return GetBodyPosVel(t_tdb, frame_centers.at(frame), target, frame, units, scale);
  }

  Vec3 GetBodyPos(Real t_tdb, BodyId target, Frame frame) {
    LUPNT_CHECK(frame != Frame::UNDEFINED, "Frame is undefined", "GetBodyPos");
    return GetBodyPos(t_tdb, frame_centers.at(frame), target, frame);
  }

  Vec3 GetBodyPos(Real t_tdb, BodyId center, BodyId target, Frame frame, const UnitSystem& units,
                  CoordinateScale scale) {
    Vec3 r = GetBodyPos(t_tdb, center, target, frame);
    r = ScalePositionForCoordinateScale(r, CoordinateScale::TDB, scale);
    return PositionFromSI(r, units);
  }

  Vec3 GetBodyPos(Real t_tdb, BodyId target, Frame frame, const UnitSystem& units,
                  CoordinateScale scale) {
    LUPNT_CHECK(frame != Frame::UNDEFINED, "Frame is undefined", "GetBodyPos");
    return GetBodyPos(t_tdb, frame_centers.at(frame), target, frame, units, scale);
  }

  MatX6 GetBodyPosVel(const VecX& t_tdb, BodyId target, Frame frame) {
    MatX6 rv(t_tdb.size(), 6);
    for (int i = 0; i < t_tdb.size(); i++) {
      rv.row(i) = GetBodyPosVel(t_tdb(i), target, frame);
    }
    return rv;
  }

  MatX3 GetBodyPos(const VecX& t_tdb, BodyId target, Frame frame) {
    MatX3 r(t_tdb.size(), 3);
    for (int i = 0; i < t_tdb.size(); i++) {
      r.row(i) = GetBodyPos(t_tdb(i), target, frame);
    }
    return r;
  }

  MatX6 GetBodyPosVel(const VecX& t_tdb, BodyId center, BodyId target, Frame frame) {
    MatX6 rv(t_tdb.size(), 6);
    for (int i = 0; i < t_tdb.size(); i++) {
      rv.row(i) = GetBodyPosVel(t_tdb(i), center, target, frame);
    }
    return rv;
  }

  MatX6 GetBodyPosVel(const VecX& t_tdb, BodyId target, Frame frame, const UnitSystem& units,
                      CoordinateScale scale) {
    MatX6 rv(t_tdb.size(), 6);
    for (int i = 0; i < t_tdb.size(); i++) {
      rv.row(i) = GetBodyPosVel(t_tdb(i), target, frame, units, scale);
    }
    return rv;
  }

  MatX6 GetBodyPosVel(const VecX& t_tdb, BodyId center, BodyId target, Frame frame,
                      const UnitSystem& units, CoordinateScale scale) {
    MatX6 rv(t_tdb.size(), 6);
    for (int i = 0; i < t_tdb.size(); i++) {
      rv.row(i) = GetBodyPosVel(t_tdb(i), center, target, frame, units, scale);
    }
    return rv;
  }

  MatX3 GetBodyPos(const VecX& t_tdb, BodyId center, BodyId target, Frame frame) {
    MatX3 r(t_tdb.size(), 3);
    for (int i = 0; i < t_tdb.size(); i++) {
      r.row(i) = GetBodyPos(t_tdb(i), center, target, frame);
    }
    return r;
  }

  MatX3 GetBodyPos(const VecX& t_tdb, BodyId target, Frame frame, const UnitSystem& units,
                   CoordinateScale scale) {
    MatX3 r(t_tdb.size(), 3);
    for (int i = 0; i < t_tdb.size(); i++) {
      r.row(i) = GetBodyPos(t_tdb(i), target, frame, units, scale);
    }
    return r;
  }

  MatX3 GetBodyPos(const VecX& t_tdb, BodyId center, BodyId target, Frame frame,
                   const UnitSystem& units, CoordinateScale scale) {
    MatX3 r(t_tdb.size(), 3);
    for (int i = 0; i < t_tdb.size(); i++) {
      r.row(i) = GetBodyPos(t_tdb(i), center, target, frame, units, scale);
    }
    return r;
  }

  double GetTtTdbDifference(double t_tdb) {
    (void)t_tdb;
    throw std::runtime_error("Not implemented");
  }

}  // namespace lupnt
