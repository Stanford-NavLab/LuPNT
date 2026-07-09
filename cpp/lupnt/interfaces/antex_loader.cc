/**
 * @file antex_loader.cc
 * @author Stanford NAV LAB
 * @brief ANTEX (antenna exchange format) file loader
 * @version 0.1
 * @date 2025-06-07
 *
 * @copyright Copyright (c) 2025
 */
#include "lupnt/interfaces/antex_loader.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <regex>
#include <sstream>

#include "lupnt/conversions/time_conversions.h"
#include "lupnt/core/error.h"
#include "lupnt/core/file.h"
#include "lupnt/interfaces/kernels.h"

namespace lupnt {

  namespace {
    constexpr int kLabelCol = 60;  // ANTEX label starts at column 61 (1-based)

    std::string GetLabel(const std::string& line) {
      if (line.size() <= kLabelCol) return "";
      std::string label = line.substr(kLabelCol);
      // rstrip
      label.erase(label.find_last_not_of(" \t\r\n") + 1);
      return label;
    }

    std::string GetData(const std::string& line) {
      return line.substr(0, std::min(line.size(), static_cast<size_t>(kLabelCol)));
    }

    std::string ToUpper(const std::string& s) {
      std::string out = s;
      std::transform(out.begin(), out.end(), out.begin(),
                     [](unsigned char c) { return std::toupper(c); });
      return out;
    }

    std::string Trim(const std::string& s) {
      size_t b = s.find_first_not_of(" \t\r\n");
      if (b == std::string::npos) return "";
      size_t e = s.find_last_not_of(" \t\r\n");
      return s.substr(b, e - b + 1);
    }

    /// @brief Parse "yyyy mm dd hh mm ss.sss ... VALID FROM/UNTIL" data field
    /// into TAI seconds (treating the timestamp as GPS time, matching the
    /// convention used elsewhere in these loaders -- the exact time-system
    /// identity is immaterial given ANTEX validity spans are weeks-to-years
    /// wide).
    bool ParseAntexDatetimeTai(const std::string& data, double& t_tai) {
      std::istringstream ss(data);
      int y, mo, d;
      int hh = 0, mm = 0;
      double ssec = 0.0;
      if (!(ss >> y >> mo >> d)) return false;
      ss >> hh >> mm >> ssec;  // optional fields default to 0

      Real t = ConvertTime(GregorianToTime(y, mo, d, hh, mm, ssec), Time::GPS, Time::TAI);
      t_tai = t.val();
      return true;
    }

    bool LooksLikePrnToken(const std::string& tok) {
      static const std::regex re("^[A-Z][0-9]{2}$");
      return std::regex_match(tok, re);
    }

    inline Vec3d ToVec3d(const Vec3& v) { return Vec3d(v(0).val(), v(1).val(), v(2).val()); }
  }  // namespace

  // Constructors ***************************************************************

  AntexLoader::AntexLoader(const std::filesystem::path& filepath) { LoadFile(filepath); }

  // Loading ********************************************************************

  void AntexLoader::LoadFile(const std::filesystem::path& filepath) { ParseFile(filepath); }

  void AntexLoader::ParseFile(const std::filesystem::path& filepath) {
    LUPNT_CHECK(std::filesystem::exists(filepath), "ANTEX file not found: " + filepath.string(),
                "AntexLoader");

    std::ifstream file = OpenFile<std::ifstream>(filepath);

    SatAntennaEntry cur_entry;
    bool have_entry = false;
    std::string cur_freq;
    bool have_freq = false;

    std::string raw;
    while (std::getline(file, raw)) {
      // strip trailing CR (Windows line endings)
      if (!raw.empty() && raw.back() == '\r') raw.pop_back();

      std::string label = GetLabel(raw);
      std::string data = GetData(raw);

      if (label == "START OF ANTENNA") {
        cur_entry = SatAntennaEntry{};
        have_entry = true;
        cur_freq.clear();
        have_freq = false;
        continue;
      }

      if (label == "END OF ANTENNA") {
        if (have_entry && !cur_entry.sat_id.empty()) {
          sat_index_[cur_entry.sat_id].push_back(cur_entry);
        }
        have_entry = false;
        cur_freq.clear();
        have_freq = false;
        continue;
      }

      if (!have_entry) continue;

      if (label == "TYPE / SERIAL NO") {
        std::istringstream ss(data);
        std::string tok;
        std::string prn_token;
        while (ss >> tok) {
          if (LooksLikePrnToken(tok)) {
            prn_token = tok;
            break;
          }
        }
        cur_entry.sat_id = prn_token;  // empty -> receiver block, not indexed
        continue;
      }

      if (label == "VALID FROM") {
        double t;
        if (ParseAntexDatetimeTai(data, t)) {
          cur_entry.has_valid_from = true;
          cur_entry.valid_from_tai = t;
        }
        continue;
      }

      if (label == "VALID UNTIL") {
        double t;
        if (ParseAntexDatetimeTai(data, t)) {
          cur_entry.has_valid_until = true;
          cur_entry.valid_until_tai = t;
        }
        continue;
      }

      if (label == "START OF FREQUENCY") {
        if (cur_entry.sat_id.empty()) {
          have_freq = false;
          continue;
        }
        std::istringstream ss(data);
        std::string tok;
        if (!(ss >> tok)) continue;
        cur_freq = Trim(tok);
        have_freq = true;
        cur_entry.freqs.emplace(cur_freq, FreqPattern{cur_freq});
        continue;
      }

      if (label == "END OF FREQUENCY") {
        cur_freq.clear();
        have_freq = false;
        continue;
      }

      if (cur_entry.sat_id.empty() || !have_freq) continue;

      if (label == "NORTH / EAST / UP") {
        std::istringstream ss(data);
        double n, e, u;
        if (ss >> n >> e >> u) {
          FreqPattern& pat = cur_entry.freqs.at(cur_freq);
          pat.has_pco = true;
          // ANTEX stores PCO in millimeters
          pat.pco_neu_m = Vec3d(n, e, u) / 1000.0;
        }
        continue;
      }
    }
  }

  // Queries ********************************************************************

  std::string AntexLoader::GnssLetter(GnssConst gnss_const) {
    switch (gnss_const) {
      case GnssConst::GPS: return "G";
      case GnssConst::GLONASS: return "R";
      case GnssConst::GALILEO: return "E";
      case GnssConst::BEIDOU: return "C";
      case GnssConst::QZSS: return "J";
      default: return "G";
    }
  }

  std::string AntexLoader::SatId(GnssConst gnss_const, int prn) {
    std::ostringstream ss;
    ss << GnssLetter(gnss_const) << std::setfill('0') << std::setw(2) << prn;
    return ss.str();
  }

  std::string AntexLoader::FreqToAntexCode(const std::string& gnss_letter,
                                           const std::string& freq_name) {
    std::string c = ToUpper(gnss_letter);
    std::string f = ToUpper(freq_name);

    if (c == "G") {
      if (f == "L1") return "G01";
      if (f == "L2") return "G02";
      if (f == "L5") return "G05";
    } else if (c == "E") {
      if (f == "E1") return "E01";
      if (f == "E5A" || f == "E5a") return "E05";
      if (f == "E5B") return "E07";
      if (f == "E5") return "E08";
      if (f == "E6") return "E06";
    } else if (c == "R") {
      if (f == "G1") return "R01";
      if (f == "G2") return "R02";
      if (f == "G3") return "R03";
    } else if (c == "J") {
      if (f == "L1") return "J01";
      if (f == "L2") return "J02";
      if (f == "L5") return "J05";
      if (f == "L6") return "J06";
    } else if (c == "C") {
      if (f == "B1I") return "C02";
      if (f == "B1C") return "C01";
      if (f == "B2A") return "C05";
      if (f == "B2I") return "C07";
      if (f == "B3I") return "C06";
    }
    return f;  // fallback: assume `freq_name` is already an ANTEX code
  }

  bool AntexLoader::HasSatellite(GnssConst gnss_const, int prn) const {
    return sat_index_.find(SatId(gnss_const, prn)) != sat_index_.end();
  }

  const AntexLoader::SatAntennaEntry& AntexLoader::SelectEntry(const std::string& sat_id,
                                                               double t_tai) const {
    auto it = sat_index_.find(sat_id);
    LUPNT_CHECK(it != sat_index_.end(), "Satellite '" + sat_id + "' not found in ANTEX file",
                "AntexLoader");

    const SatAntennaEntry* best = nullptr;
    for (const auto& e : it->second) {
      bool ok_from = !e.has_valid_from || t_tai >= e.valid_from_tai;
      bool ok_until = !e.has_valid_until || t_tai <= e.valid_until_tai;
      if (ok_from && ok_until) {
        // Prefer the entry with the latest `valid_from` (mirrors `_select_entry`'s
        // descending sort by `valid_from`).
        if (best == nullptr
            || (e.has_valid_from
                && (!best->has_valid_from || e.valid_from_tai > best->valid_from_tai))) {
          best = &e;
        }
      }
    }
    LUPNT_CHECK(best != nullptr,
                "No valid ANTEX entry for satellite '" + sat_id + "' at the requested epoch",
                "AntexLoader");
    return *best;
  }

  std::vector<std::string> AntexLoader::GetAvailableFreqCodes(GnssConst gnss_const, int prn,
                                                              Real t_tai) const {
    const SatAntennaEntry& entry = SelectEntry(SatId(gnss_const, prn), t_tai.val());
    std::vector<std::string> codes;
    codes.reserve(entry.freqs.size());
    for (const auto& [code, _] : entry.freqs) codes.push_back(code);
    std::sort(codes.begin(), codes.end());
    return codes;
  }

  Vec3d AntexLoader::GetPco(const std::string& gnss_letter, int prn,
                            const std::string& antex_freq_code, Real t_tai) const {
    std::ostringstream ss;
    ss << ToUpper(gnss_letter) << std::setfill('0') << std::setw(2) << prn;
    std::string sat_id = ss.str();

    const SatAntennaEntry& entry = SelectEntry(sat_id, t_tai.val());
    std::string code = ToUpper(antex_freq_code);
    auto it = entry.freqs.find(code);
    LUPNT_CHECK(it != entry.freqs.end() && it->second.has_pco,
                "PCO not found for satellite '" + sat_id + "', frequency code '" + code + "'",
                "AntexLoader");
    return it->second.pco_neu_m;
  }

  Vec3d AntexLoader::GetPco(GnssConst gnss_const, int prn, GnssFreq freq, Real t_tai) const {
    std::string gnss_letter = GnssLetter(gnss_const);
    std::string freq_name = ToUpper(std::string(enum_name(freq)));
    std::string code = FreqToAntexCode(gnss_letter, freq_name);
    return GetPco(gnss_letter, prn, code, t_tai);
  }

  bool AntexLoader::HasPco(GnssConst gnss_const, int prn, GnssFreq freq, Real t_tai) const {
    if (!HasSatellite(gnss_const, prn)) return false;
    std::string gnss_letter = GnssLetter(gnss_const);
    std::string freq_name = ToUpper(std::string(enum_name(freq)));
    std::string code = ToUpper(FreqToAntexCode(gnss_letter, freq_name));
    std::ostringstream ss;
    ss << ToUpper(gnss_letter) << std::setfill('0') << std::setw(2) << prn;
    const SatAntennaEntry& entry = SelectEntry(ss.str(), t_tai.val());
    auto it = entry.freqs.find(code);
    return it != entry.freqs.end() && it->second.has_pco;
  }

  // PCO -> ECEF correction *****************************************************

  Mat3d AntexLoader::ComputeIjkToEcefRotation(Real t_tai, const Vec3d& r_sat_ecef) {
    constexpr double kEps = 1e-12;

    Vec3d kvec = -r_sat_ecef.normalized();
    LUPNT_CHECK(r_sat_ecef.norm() > kEps, "Cannot normalize near-zero satellite position",
                "AntexLoader");

    Real t_tdb = ConvertTime(t_tai, Time::TAI, Time::TDB);
    Vec3d r_sun_ecef = ToVec3d(GetBodyPos(t_tdb, BodyId::EARTH, BodyId::SUN, Frame::ECEF));
    Vec3d to_sun = r_sun_ecef - r_sat_ecef;
    LUPNT_CHECK(to_sun.norm() > kEps, "Cannot normalize near-zero Sun direction", "AntexLoader");
    Vec3d jvec = to_sun.normalized();

    Vec3d ivec = jvec.cross(kvec);

    Mat3d Cijk;
    Cijk.col(0) = ivec;
    Cijk.col(1) = jvec;
    Cijk.col(2) = kvec;
    return Cijk;
  }

  Vec3d AntexLoader::ApplyPcoCorrectionEcef(Real t_tai, const Vec3d& pos_sp3_ecef,
                                            const Vec3d& pco_neu_m) {
    Mat3d Cijk = ComputeIjkToEcefRotation(t_tai, pos_sp3_ecef);
    Vec3d pco_ecef = Cijk * pco_neu_m;
    return pos_sp3_ecef + pco_ecef;
  }

}  // namespace lupnt
