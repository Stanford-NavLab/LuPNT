#include "lupnt/interfaces/eop.h"

#include <fmt/format.h>

#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "lupnt/core/constants.h"
#include "lupnt/core/error.h"
#include "lupnt/core/file.h"
#include "lupnt/core/logger.h"
#include "lupnt/numerics/interpolation.h"

namespace lupnt {

  UniquePtr<EopFileData> eop;
  std::mutex eop_mutex;

  namespace {
    // Guarded by eop_mutex. Reset on every (re)load so a fresh table gets a fresh warning.
    bool eop_clamp_warned = false;
    bool eop_prediction_warned = false;

    // The declared EOP source and optional explicit file (see SetEopSource), resolved when the
    // table is first needed. Guarded by eop_mutex. `eop_selected_path` empty means "use the
    // default file for the source"; `eop_loaded_path` records what was actually read, so
    // SetEopSource can tell whether a reload is required.
    EopSource eop_selected_source = EopSource::C04;
    EopNutation eop_selected_nutation = EopNutation::Iau1980;
    std::filesystem::path eop_selected_path;
    std::filesystem::path eop_loaded_path;

    std::filesystem::path DefaultPathFor(EopSource source) {
      switch (source) {
        case EopSource::Finals: return GetDataPath() / "planetary_coeff" / EOP_FINALS_FILENAME;
        case EopSource::C04:
        default: return GetFilePath(EOP_FILENAME);
      }
    }

    // EOP perturbation state (see SetEopPerturbation). Guarded by its own mutex, always acquired
    // *inside* eop_mutex when both are held, so the order is consistent and cannot deadlock.
    std::mutex eop_perturb_mutex;
    EopPerturbation eop_perturbation;
    std::function<EopPerturbation(Real)> eop_perturbation_fn;
    bool eop_perturbation_set = false;

    // True when the celestial-pole offsets can be non-zero -- an IAU2000 table is loaded, or a
    // perturbation is active. Atomic and checked without a lock so RotPrecessionNutation can skip
    // an EOP lookup entirely on the default (C04, unperturbed) path.
    std::atomic<bool> eop_cpo_active{false};

    bool CpoFromTableLocked() { return eop && eop->nutation == EopNutation::Iau2000; }

    // Last epoch in the loaded table backed by a final/measured value. Cached at load time rather
    // than recomputed per call: GetEopData sits in the frame-conversion inner loop and the scan
    // is O(rows) over ~20k rows. Guarded by eop_mutex.
    double eop_last_measured_mjd = 0.0;

    /// Scan the loaded table for the last non-predicted row. Caller must hold eop_mutex with
    /// `eop` already populated; call once per load and cache in eop_last_measured_mjd.
    double ScanLastMeasuredMjdLocked() {
      Eigen::Index n = eop->mjds_utc.size();
      if (eop->is_prediction.size() != n) return eop->mjds_utc(n - 1);
      // Start below the table so an all-predicted table reports "nothing measured" rather than
      // silently claiming the first (or last) row as final.
      double last = eop->mjds_utc(0) - 1.0;
      for (Eigen::Index i = 0; i < n; ++i) {
        if (eop->is_prediction(i) == 0) last = eop->mjds_utc(i);
      }
      return last;
    }

    /// Download `url` to `dest_path` using curl. The file is first written to a
    /// temporary path and atomically moved into place on success, leaving any
    /// existing file at `dest_path` untouched on failure.
    bool DownloadFile(const std::string& url, const std::filesystem::path& dest_path) {
      std::filesystem::path tmp_path = dest_path;
      tmp_path += ".part";
      std::error_code ec;
      std::filesystem::remove(tmp_path, ec);

      std::filesystem::create_directories(dest_path.parent_path(), ec);

      std::string cmd
          = fmt::format("curl -fsSL --connect-timeout 10 --max-time 600 -o \"{}\" \"{}\"",
                        tmp_path.string(), url);
      bool ok = std::system(cmd.c_str()) == 0 && std::filesystem::exists(tmp_path)
                && std::filesystem::file_size(tmp_path) > 0;
      if (!ok) {
        std::filesystem::remove(tmp_path, ec);
        return false;
      }
      std::filesystem::rename(tmp_path, dest_path, ec);
      if (ec) {
        ec.clear();
        std::filesystem::copy_file(tmp_path, dest_path,
                                   std::filesystem::copy_options::overwrite_existing, ec);
        std::filesystem::remove(tmp_path, ec);
      }
      return !ec;
    }
    /// Warn when a load is dropped because a table is already present. This is the one way to end
    /// up silently running on a different EOP product than the one just asked for, so it must not
    /// pass unremarked; SetEopSource is the ordering-independent way to express the choice.
    void WarnIgnoredLoadLocked(const std::filesystem::path& filepath) {
      Logger::Warn(
          fmt::format("EOP: ignoring load of '{}' -- a {} table from '{}' is already loaded. Pass "
                      "force=true to replace it, or use SetEopSource() before first use.",
                      filepath.filename().string(),
                      eop->source == EopSource::Finals ? "Finals" : "C04",
                      eop_loaded_path.filename().string()),
          "Eop");
    }

    // Unlocked core of LoadEopFileData -- caller must already hold eop_mutex. Split out so
    // GetEopFileData/GetEopData can perform their "is it loaded yet" check and the load itself
    // under a single critical section (see the comment on those functions for why the previous
    // unlocked check-then-load-then-read was a data race).
    void LoadEopFileDataLocked(const std::filesystem::path& filepath, bool force) {
      if (eop && !force) {  // Data already loaded
        WarnIgnoredLoadLocked(filepath);
        return;
      }

      int n_header_lines = 14;
      size_t n_lines = CountLines(filepath.string()) - n_header_lines;
      std::ifstream file = OpenFile<std::ifstream>(filepath);
      // Skip header lines
      std::string line;
      for (int i = 0; i < n_header_lines; ++i) {
        std::getline(file, line);
      }

      // Initialize EopFileData struct
      eop = MakeUnique<EopFileData>();
      eop_clamp_warned = false;
      eop_prediction_warned = false;
      eop->years.resize(n_lines);
      eop->months.resize(n_lines);
      eop->days.resize(n_lines);
      eop->mjds_utc.resize(n_lines);
      eop->x.resize(n_lines);
      eop->y.resize(n_lines);
      eop->ut1_utc.resize(n_lines);
      eop->lod.resize(n_lines);
      eop->dpsi.resize(n_lines);
      eop->deps.resize(n_lines);
      eop->xErr.resize(n_lines);
      eop->yErr.resize(n_lines);
      eop->ut1_utc_err.resize(n_lines);
      eop->lod_err.resize(n_lines);
      eop->dpsi_err.resize(n_lines);
      eop->deps_err.resize(n_lines);
      eop->is_prediction.resize(n_lines);
      eop->source = EopSource::C04;
      // C04 (IAU1980) carries dpsi/deps, not CIP offsets. Zero-fill dX/dY rather than leaving
      // them empty so every consumer can index them unconditionally.
      eop->nutation = EopNutation::Iau1980;
      eop->dX = VecXd::Zero(n_lines);
      eop->dY = VecXd::Zero(n_lines);
      eop->dX_err = VecXd::Zero(n_lines);
      eop->dY_err = VecXd::Zero(n_lines);

      size_t row = 0;
      // Read data lines
      while (std::getline(file, line)) {
        if (line.empty()) {
          continue;
        }
        std::istringstream iss(line);
        int year, month, day, mjd;
        double x, y, ut1_utc, lod, dpsi, deps, xErr, yErr, ut1_utc_err, lod_err, dpsi_err, deps_err;

        iss >> year >> month >> day >> mjd >> x >> y >> ut1_utc >> lod >> dpsi >> deps >> xErr
            >> yErr >> ut1_utc_err >> lod_err >> dpsi_err >> deps_err;

        eop->years(row) = year;
        eop->months(row) = month;
        eop->days(row) = day;
        eop->mjds_utc(row) = mjd;
        eop->x(row) = x;
        eop->y(row) = y;
        eop->ut1_utc(row) = ut1_utc;
        eop->lod(row) = lod;
        eop->dpsi(row) = dpsi;
        eop->deps(row) = deps;
        eop->xErr(row) = xErr;
        eop->yErr(row) = yErr;
        eop->ut1_utc_err(row) = ut1_utc_err;
        eop->lod_err(row) = lod_err;
        eop->dpsi_err(row) = dpsi_err;
        eop->deps_err(row) = deps_err;
        // C04 marks rows whose values are not final by filling every formal-error column with
        // 0.999. Those rows are not merely missing their sigmas -- their EOP values are
        // preliminary/extrapolated and can differ from the settled series by tens of mas in
        // polar motion, so they must not be treated as measured.
        eop->is_prediction(row) = xErr >= 0.9 ? 1 : 0;

        ++row;
      }

      file.close();
      eop_last_measured_mjd = ScanLastMeasuredMjdLocked();
      eop_selected_source = EopSource::C04;
      eop_selected_nutation = EopNutation::Iau1980;
      eop_loaded_path = filepath;
      eop_cpo_active = CpoFromTableLocked() || eop_perturbation_set;
      return;
    }

    /// Extract a fixed-width field from an IERS finals line. Returns nullopt when the field is
    /// blank or the line is too short: finals rows drop trailing columns as data runs out, so a
    /// blank is normal and must not be read as a zero.
    std::optional<double> ParseFinalsField(const std::string& line, size_t begin, size_t end) {
      if (line.size() <= begin) return std::nullopt;
      std::string field = line.substr(begin, std::min(end, line.size()) - begin);
      size_t first = field.find_first_not_of(" \t\r\n");
      if (first == std::string::npos) return std::nullopt;
      size_t last = field.find_last_not_of(" \t\r\n");
      try {
        return std::stod(field.substr(first, last - first + 1));
      } catch (const std::exception&) {
        return std::nullopt;
      }
    }

    char FinalsFlag(const std::string& line, size_t idx) {
      return line.size() > idx ? line[idx] : ' ';
    }

    // Unlocked core of LoadEopFinalsFileData -- caller must already hold eop_mutex, for the same
    // reason as LoadEopFileDataLocked.
    //
    // Column positions are the documented IERS finals format (1-indexed in the IERS readme,
    // 0-indexed half-open here):
    //   yy 0-2, mm 2-4, dd 4-6, MJD 7-15, PM flag 16, x 18-27, sig_x 27-36, y 37-46,
    //   sig_y 46-55, UT1 flag 57, UT1-UTC 58-68, sig_UT1 68-78, LOD 79-86, sig_LOD 86-93,
    //   nutation flag 95, dPsi 96-106, sig_dPsi 106-116, dEps 116-125, sig_dEps 125-134.
    // Units differ from C04: LOD and the nutation columns are in milliseconds and milliarcseconds
    // respectively, and are scaled here to seconds and arcseconds to match EopFileData.
    void LoadEopFinalsFileDataLocked(const std::filesystem::path& filepath, bool force,
                                     EopNutation nutation) {
      if (eop && !force) {  // Data already loaded
        WarnIgnoredLoadLocked(filepath);
        return;
      }

      std::ifstream file = OpenFile<std::ifstream>(filepath);
      std::string line;

      std::vector<int> years, months, days, is_prediction;
      std::vector<double> mjds, x, y, ut1_utc, lod, dpsi, deps;
      std::vector<double> x_err, y_err, ut1_utc_err, lod_err, dpsi_err, deps_err;

      // Optional columns are held at their last valid value once the file stops carrying them,
      // and where each ran out is reported below so the hold is visible rather than assumed.
      double held_lod = 0.0, held_lod_err = 0.0;
      double held_dpsi = 0.0, held_deps = 0.0, held_dpsi_err = 0.0, held_deps_err = 0.0;
      double mjd_lod_end = 0.0, mjd_nutation_end = 0.0;
      bool have_lod = false, have_nutation = false;

      while (std::getline(file, line)) {
        std::optional<double> mjd = ParseFinalsField(line, 7, 15);
        std::optional<double> xv = ParseFinalsField(line, 18, 27);
        std::optional<double> yv = ParseFinalsField(line, 37, 46);
        std::optional<double> uv = ParseFinalsField(line, 58, 68);
        // Polar motion and UT1-UTC are the essential trio; without them the row carries no usable
        // orientation (the file ends with a run of date-only rows past the prediction horizon).
        if (!mjd || !xv || !yv || !uv) continue;

        std::optional<double> lodv = ParseFinalsField(line, 79, 86);
        if (lodv) {
          held_lod = *lodv * 1e-3;  // ms -> s
          held_lod_err = ParseFinalsField(line, 86, 93).value_or(0.0) * 1e-3;
          mjd_lod_end = *mjd;
          have_lod = true;
        }
        std::optional<double> dpsiv = ParseFinalsField(line, 96, 106);
        std::optional<double> depsv = ParseFinalsField(line, 116, 125);
        if (dpsiv && depsv) {
          held_dpsi = *dpsiv * 1e-3;  // mas -> arcsec
          held_deps = *depsv * 1e-3;
          held_dpsi_err = ParseFinalsField(line, 106, 116).value_or(0.0) * 1e-3;
          held_deps_err = ParseFinalsField(line, 125, 134).value_or(0.0) * 1e-3;
          mjd_nutation_end = *mjd;
          have_nutation = true;
        }

        // Two-digit year; the finals series starts in 1973 and the window is unambiguous well
        // past this code's lifetime.
        int yy = static_cast<int>(ParseFinalsField(line, 0, 2).value_or(0.0));
        years.push_back(yy >= 73 ? 1900 + yy : 2000 + yy);
        months.push_back(static_cast<int>(ParseFinalsField(line, 2, 4).value_or(0.0)));
        days.push_back(static_cast<int>(ParseFinalsField(line, 4, 6).value_or(0.0)));

        mjds.push_back(*mjd);
        x.push_back(*xv);
        y.push_back(*yv);
        ut1_utc.push_back(*uv);
        lod.push_back(held_lod);
        dpsi.push_back(held_dpsi);
        deps.push_back(held_deps);
        x_err.push_back(ParseFinalsField(line, 27, 36).value_or(0.0));
        y_err.push_back(ParseFinalsField(line, 46, 55).value_or(0.0));
        ut1_utc_err.push_back(ParseFinalsField(line, 68, 78).value_or(0.0));
        lod_err.push_back(held_lod_err);
        dpsi_err.push_back(held_dpsi_err);
        deps_err.push_back(held_deps_err);

        // The polar-motion and UT1 flags are independent; a row is predicted if either is.
        // (The nutation flag is deliberately excluded -- it reads `P` even on rows whose polar
        // motion and UT1 are final.)
        bool predicted = FinalsFlag(line, 16) == 'P' || FinalsFlag(line, 57) == 'P';
        is_prediction.push_back(predicted ? 1 : 0);
      }
      file.close();

      LUPNT_CHECK(!mjds.empty(),
                  fmt::format("EOP: no usable data rows found in finals file '{}' -- is it "
                              "actually an IERS finals file?",
                              filepath.string()),
                  "Eop");

      auto to_vecxd = [](const std::vector<double>& v) {
        return VecXd(Eigen::Map<const VecXd>(v.data(), static_cast<Eigen::Index>(v.size())));
      };
      auto to_vecxi = [](const std::vector<int>& v) {
        return VecXi(Eigen::Map<const VecXi>(v.data(), static_cast<Eigen::Index>(v.size())));
      };

      eop = MakeUnique<EopFileData>();
      eop_clamp_warned = false;
      eop_prediction_warned = false;
      eop->source = EopSource::Finals;
      eop->years = to_vecxi(years);
      eop->months = to_vecxi(months);
      eop->days = to_vecxi(days);
      eop->mjds_utc = to_vecxd(mjds);
      eop->x = to_vecxd(x);
      eop->y = to_vecxd(y);
      eop->ut1_utc = to_vecxd(ut1_utc);
      eop->lod = to_vecxd(lod);
      eop->dpsi = to_vecxd(dpsi);
      eop->deps = to_vecxd(deps);
      eop->xErr = to_vecxd(x_err);
      eop->yErr = to_vecxd(y_err);
      eop->ut1_utc_err = to_vecxd(ut1_utc_err);
      eop->lod_err = to_vecxd(lod_err);
      // The two finals variants share a layout but not a meaning: the same columns are dPsi/dEps
      // in finals.all and the CIP offsets dX/dY in finals2000A.all. File them under whichever the
      // caller declared and zero-fill the other pair.
      eop->nutation = nutation;
      VecXd zeros = VecXd::Zero(static_cast<Eigen::Index>(mjds.size()));
      if (nutation == EopNutation::Iau2000) {
        eop->dX = to_vecxd(dpsi);
        eop->dY = to_vecxd(deps);
        eop->dX_err = to_vecxd(dpsi_err);
        eop->dY_err = to_vecxd(deps_err);
        eop->dpsi = zeros;
        eop->deps = zeros;
        eop->dpsi_err = zeros;
        eop->deps_err = zeros;
      } else {
        eop->dpsi_err = to_vecxd(dpsi_err);
        eop->deps_err = to_vecxd(deps_err);
        eop->dX = zeros;
        eop->dY = zeros;
        eop->dX_err = zeros;
        eop->dY_err = zeros;
      }
      eop->is_prediction = to_vecxi(is_prediction);

      eop_last_measured_mjd = ScanLastMeasuredMjdLocked();
      eop_selected_source = EopSource::Finals;
      eop_selected_nutation = nutation;
      eop_loaded_path = filepath;
      eop_cpo_active = CpoFromTableLocked() || eop_perturbation_set;
      double mjd_last_measured = eop_last_measured_mjd;
      Logger::Info(
          fmt::format("EOP: loaded {} finals rows from '{}', MJD [{:.1f}, {:.1f}]; measured "
                      "through {:.1f}, predicted for {:.0f} days beyond that.",
                      mjds.size(), filepath.filename().string(), mjds.front(), mjds.back(),
                      mjd_last_measured, mjds.back() - mjd_last_measured),
          "Eop");
      if (have_lod && mjd_lod_end < mjds.back()) {
        Logger::Info(fmt::format("EOP: finals LOD ends at MJD {:.1f} and the nutation columns at "
                                 "MJD {:.1f}; the last valid value of each is held constant "
                                 "beyond those epochs.",
                                 mjd_lod_end, have_nutation ? mjd_nutation_end : mjd_lod_end),
                     "Eop");
      }
    }

    /// Add the active perturbation to an interpolated EopData. Called for every GetEopData
    /// result, so a perturbation reaches polar motion, the UT1 time conversion (and hence
    /// sidereal rotation), LOD and the CIP offsets through one seam instead of a parallel stack.
    void ApplyEopPerturbationLocked(Real mjd_utc, EopData& data) {
      if (!eop_perturbation_set) return;
      std::lock_guard<std::mutex> lock(eop_perturb_mutex);
      EopPerturbation p = eop_perturbation_fn ? eop_perturbation_fn(mjd_utc) : eop_perturbation;
      data.x_pole += p.dx_pole;
      data.y_pole += p.dy_pole;
      data.ut1_utc += p.dut1;
      data.lod += p.dlod;
      data.dX += p.ddX;
      data.dY += p.ddY;
    }

    /// Load the declared EOP source if nothing is loaded yet. This is the single place the
    /// lazy-load default is decided; every entry point that needs the table routes through it, so
    /// the source cannot depend on which caller happened to arrive first.
    void EnsureEopLoadedLocked() {
      if (eop) return;
      std::filesystem::path path
          = eop_selected_path.empty() ? DefaultPathFor(eop_selected_source) : eop_selected_path;
      if (eop_selected_source == EopSource::Finals) {
        LoadEopFinalsFileDataLocked(path, false, eop_selected_nutation);
      } else {
        LoadEopFileDataLocked(path, false);
      }
    }

  }  // namespace

  void SetEopSource(EopSource source, const std::filesystem::path& filepath) {
    std::lock_guard<std::mutex> lock(eop_mutex);
    std::filesystem::path path = filepath.empty() ? DefaultPathFor(source) : filepath;

    // Fail here, at the point of the mistake, rather than deep inside the first frame conversion
    // that happens to need Earth orientation. Nothing is bundled for Finals, so a missing file is
    // the expected first-run outcome and the message has to say what to do about it.
    LUPNT_CHECK(std::filesystem::exists(path),
                fmt::format("EOP: cannot select {} source -- file '{}' does not exist.{}",
                            source == EopSource::Finals ? "Finals" : "C04", path.string(),
                            source == EopSource::Finals
                                ? " Call LoadLatestEopFinalsFromIers() to download one, or pass"
                                  " an explicit path."
                                : ""),
                "Eop");

    eop_selected_source = source;
    eop_selected_path = filepath;  // stays empty when defaulted, so the default can move

    // A table already being loaded is exactly the case this API exists for: deferring here would
    // silently ignore the call, which is the failure mode SetEopSource is meant to remove.
    if (eop && (eop->source != source || eop_loaded_path != path)) {
      Logger::Info(
          fmt::format("EOP: switching source to {} ('{}'), replacing the loaded table.",
                      source == EopSource::Finals ? "Finals" : "C04", path.filename().string()),
          "Eop");
      if (source == EopSource::Finals) {
        LoadEopFinalsFileDataLocked(path, true, eop_selected_nutation);
      } else {
        LoadEopFileDataLocked(path, true);
      }
    }
  }

  EopSource GetEopSource() {
    std::lock_guard<std::mutex> lock(eop_mutex);
    return eop_selected_source;
  }

  void LoadEopFileData(const std::filesystem::path& filepath, bool force) {
    std::lock_guard<std::mutex> lock(eop_mutex);
    LoadEopFileDataLocked(filepath, force);
    // Only when the load actually happened -- with force=false it may have been a no-op.
    if (eop_loaded_path == filepath) eop_selected_path = filepath;
  }

  bool LoadLatestEopFromIers(bool force) {
    static const std::string url
        = "https://datacenter.iers.org/data/latestVersion/EOP_14_C04_IAU1980_one_file_1962-now.txt";
    std::filesystem::path dest_path = GetDataPath() / "planetary_coeff" / "eopc04_iers_latest.txt";

    if (DownloadFile(url, dest_path)) {
      Logger::Info(fmt::format("EOP: downloaded latest IERS EOP data to '{}'", dest_path.string()),
                   "Eop");
      LoadEopFileData(dest_path, force);
      return true;
    }

    Logger::Warn("EOP: failed to download latest IERS EOP data, falling back to bundled data",
                 "Eop");
    GetEopFileData();  // Ensure the bundled file is loaded as a fallback
    return false;
  }

  void LoadEopFinalsFileData(const std::filesystem::path& filepath, bool force,
                             EopNutation nutation) {
    std::lock_guard<std::mutex> lock(eop_mutex);
    LoadEopFinalsFileDataLocked(filepath, force, nutation);
    if (eop_loaded_path == filepath) eop_selected_path = filepath;
  }

  bool LoadLatestEopFinalsFromIers(bool force) {
    static const std::string url
        = "https://datacenter.iers.org/data/latestVersion/finals.all.iau1980.txt";
    std::filesystem::path dest_path = GetDataPath() / "planetary_coeff" / "finals.all.iau1980.txt";

    if (DownloadFile(url, dest_path)) {
      Logger::Info(
          fmt::format("EOP: downloaded latest IERS finals data to '{}'", dest_path.string()),
          "Eop");
      LoadEopFinalsFileData(dest_path, force);
      return true;
    }

    // Unlike the C04 path there is no bundled finals file to fall back to, so leave whatever is
    // loaded in place rather than silently downgrading a table the caller may have set up.
    Logger::Warn("EOP: failed to download latest IERS finals data", "Eop");
    GetEopFileData();  // Ensure some table is loaded
    return false;
  }

  Real GetUt1UtcDifference(Real mjd_utc) {
    EopData eop = GetEopData(mjd_utc);
    return eop.ut1_utc;
  }

  EopFileData* GetEopFileData() {
    // Locked for the whole check-then-load: an unlocked `if (!eop)` here raced with concurrent
    // callers (e.g. BuildConstellations's per-frequency OpenMP threads), where one thread could
    // observe a non-null `eop` pointer while another was still mid-populate inside
    // LoadEopFileDataLocked, causing reads of not-yet-resized/partially-written Eigen vectors.
    std::lock_guard<std::mutex> lock(eop_mutex);
    EnsureEopLoadedLocked();
    return eop.get();
  }

  EopData GetEopData(Real mjd_utc) {
    // Locked for the whole function (load-if-needed AND the reads below) -- see the comment on
    // GetEopFileData for why the load alone being locked isn't sufficient.
    std::lock_guard<std::mutex> lock(eop_mutex);
    EnsureEopLoadedLocked();

    EopData data;

    // Check if the requested MJD is outside the range
    if (mjd_utc <= eop->mjds_utc(0) || mjd_utc >= eop->mjds_utc(eop->mjds_utc.size() - 1)) {
      bool below = mjd_utc <= eop->mjds_utc(0);
      int i = below ? 0 : eop->mjds_utc.size() - 1;
      // Warn once per loaded table: the values below are the nearest endpoint held constant, not
      // an extrapolation, so a request far outside the table silently returns stale EOP.
      if (!eop_clamp_warned) {
        eop_clamp_warned = true;
        double gap_days = below ? eop->mjds_utc(0) - mjd_utc.val()
                                : mjd_utc.val() - eop->mjds_utc(eop->mjds_utc.size() - 1);
        Logger::Warn(
            fmt::format("EOP: epoch MJD {:.5f} is {} the loaded table [{:.1f}, {:.1f}] by {:.1f} "
                        "days; holding the nearest endpoint constant (no extrapolation). Polar "
                        "motion and UT1-UTC are stale by that much. Further clamps not reported.",
                        mjd_utc.val(), below ? "before" : "past", eop->mjds_utc(0),
                        eop->mjds_utc(eop->mjds_utc.size() - 1), gap_days),
            "Eop");
      }
      data.x_pole = eop->x(i) * RAD_ARCSEC;
      data.y_pole = eop->y(i) * RAD_ARCSEC;
      data.ut1_utc = eop->ut1_utc(i);
      data.lod = eop->lod(i);
      data.dpsi = eop->dpsi(i) * RAD_ARCSEC;
      data.deps = eop->deps(i) * RAD_ARCSEC;
      data.sigma_x_pole = eop->xErr(i) * RAD_ARCSEC;
      data.sigma_y_pole = eop->yErr(i) * RAD_ARCSEC;
      data.sigma_ut1_utc = eop->ut1_utc_err(i);
      data.dX = eop->dX(i) * RAD_ARCSEC;
      data.dY = eop->dY(i) * RAD_ARCSEC;
      ApplyEopPerturbationLocked(mjd_utc, data);
      return data;
    }

    // In range, but possibly in the predicted span. Warn once: predicted EOP carries error that
    // grows with lead time (UT1 fastest, since it integrates LOD error), which is exactly what a
    // navigation error budget needs to account for rather than inherit silently.
    double mjd_last_measured = eop_last_measured_mjd;
    if (!eop_prediction_warned && mjd_utc > mjd_last_measured) {
      eop_prediction_warned = true;
      Logger::Warn(
          fmt::format("EOP: epoch MJD {:.5f} is {:.1f} days past the last measured value "
                      "(MJD {:.1f}); using predicted EOP. Further predicted epochs not reported.",
                      mjd_utc.val(), mjd_utc.val() - mjd_last_measured, mjd_last_measured),
          "Eop");
    }

    int order = 3;
    LagrangeInterpolator interp(eop->mjds_utc, mjd_utc.val(), order);
    data.x_pole = interp.Interpolate(eop->x) * RAD_ARCSEC;
    data.y_pole = interp.Interpolate(eop->y) * RAD_ARCSEC;
    data.ut1_utc = interp.Interpolate(eop->ut1_utc);
    data.lod = interp.Interpolate(eop->lod);
    data.dpsi = interp.Interpolate(eop->dpsi) * RAD_ARCSEC;
    data.deps = interp.Interpolate(eop->deps) * RAD_ARCSEC;
    data.sigma_x_pole = interp.Interpolate(eop->xErr) * RAD_ARCSEC;
    data.sigma_y_pole = interp.Interpolate(eop->yErr) * RAD_ARCSEC;
    data.sigma_ut1_utc = interp.Interpolate(eop->ut1_utc_err);
    data.dX = interp.Interpolate(eop->dX) * RAD_ARCSEC;
    data.dY = interp.Interpolate(eop->dY) * RAD_ARCSEC;
    ApplyEopPerturbationLocked(mjd_utc, data);
    return data;
  }

  void SetEopPerturbation(const EopPerturbation& perturbation) {
    {
      std::lock_guard<std::mutex> lock(eop_perturb_mutex);
      eop_perturbation = perturbation;
      eop_perturbation_fn = nullptr;
      eop_perturbation_set = true;
    }
    eop_cpo_active = true;
  }

  void SetEopPerturbationFunction(std::function<EopPerturbation(Real)> fn) {
    {
      std::lock_guard<std::mutex> lock(eop_perturb_mutex);
      eop_perturbation_fn = std::move(fn);
      eop_perturbation = EopPerturbation{};
      eop_perturbation_set = static_cast<bool>(eop_perturbation_fn);
    }
    std::lock_guard<std::mutex> lock(eop_mutex);
    eop_cpo_active = CpoFromTableLocked() || eop_perturbation_set;
  }

  void ClearEopPerturbation() {
    {
      std::lock_guard<std::mutex> lock(eop_perturb_mutex);
      eop_perturbation = EopPerturbation{};
      eop_perturbation_fn = nullptr;
      eop_perturbation_set = false;
    }
    std::lock_guard<std::mutex> lock(eop_mutex);
    eop_cpo_active = CpoFromTableLocked();
  }

  EopPerturbation GetEopPerturbation(Real mjd_utc) {
    std::lock_guard<std::mutex> lock(eop_perturb_mutex);
    if (!eop_perturbation_set) return EopPerturbation{};
    return eop_perturbation_fn ? eop_perturbation_fn(mjd_utc) : eop_perturbation;
  }

  bool EopHasCelestialPoleOffsets() { return eop_cpo_active.load(); }

  EopCoverage GetEopCoverage() {
    std::lock_guard<std::mutex> lock(eop_mutex);
    EnsureEopLoadedLocked();
    return {eop->mjds_utc(0), eop->mjds_utc(eop->mjds_utc.size() - 1), eop_last_measured_mjd,
            eop->source};
  }

}  // namespace lupnt
