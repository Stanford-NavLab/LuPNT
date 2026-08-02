/**
 * @file SpiceInterface.cpp
 * @author Stanford NAV LAB
 * @brief  SPICE Interface functions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lupnt/interfaces/spice.h"

#include <cspice/SpiceUsr.h>
#include <cspice/SpiceZfc.h>
#include <string.h>

#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>

#include "lupnt/conversions/coordinate_conversions.h"
#include "lupnt/conversions/time_conversions.h"
#include "lupnt/core/error.h"
#include "lupnt/core/file.h"
#include "lupnt/core/logger.h"
#include "lupnt/environment/body.h"
#include "lupnt/interfaces/spice_cheby.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {

  namespace spice {
    static segment_t* cheby_s;
    static long cheby_n;

    bool spice_loaded = false;

    // True once the DE440t TT-TDB time ephemeris (NAIF body 1000000001) has
    // been loaded, enabling the high-fidelity TDB<->TT path in ConvertTime.
    static bool tt_tdb_kernel_loaded = false;

    // Base NAIF id that all time-ephemeris segments are referenced to. The
    // segment ids themselves (kNaifTtMinusTdb, ...) live in spice.h.
    static constexpr int kNaifTdb = 1000000000;

    void CheckSpiceFailure(const std::string& context) {
      if (!failed_c()) return;

      SpiceChar short_msg[1841];
      SpiceChar long_msg[1841];
      getmsg_c("SHORT", sizeof(short_msg), short_msg);
      getmsg_c("LONG", sizeof(long_msg), long_msg);
      reset_c();

      throw std::runtime_error(fmt::format("{} failed: {} {}", context, short_msg, long_msg));
    }

    void SetSpiceReturnMode() {
      SpiceChar action[] = "RETURN";
      erract_c("SET", 0, action);
    }

    namespace {
      /// Root of the NAIF generic kernels archive that hosts the Earth and
      /// lunar high-accuracy orientation kernels.
      constexpr const char* kNaifGenericKernelsUrl
          = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels";

      /// If set to a "truthy" value, disables network access during
      /// `LoadSpiceKernel` so that only locally-cached kernels are used. This
      /// is useful for offline development and for test environments without
      /// internet access.
      bool IsKernelDownloadDisabled() {
        const char* val = std::getenv("LUPNT_SKIP_SPICE_KERNEL_DOWNLOAD");
        if (val == nullptr) return false;
        std::string s(val);
        return !(s.empty() || s == "0" || s == "false" || s == "FALSE" || s == "OFF" || s == "off");
      }

      /// Run a shell command, returning true if it exits with status 0.
      bool RunShellCommand(const std::string& cmd) { return std::system(cmd.c_str()) == 0; }

      /// Download `url` to `dest_path` using curl. The file is first written
      /// to a temporary path and atomically moved into place on success, so a
      /// partially-downloaded file is never left at `dest_path`. Returns true
      /// on success.
      bool DownloadToFile(const std::string& url, const std::filesystem::path& dest_path) {
        if (IsKernelDownloadDisabled()) return false;

        std::filesystem::path tmp_path = dest_path;
        tmp_path += ".part";
        std::error_code ec;
        std::filesystem::remove(tmp_path, ec);

        std::string cmd
            = fmt::format("curl -fsSL --connect-timeout 10 --max-time 600 -o \"{}\" \"{}\"",
                          tmp_path.string(), url);
        bool ok = RunShellCommand(cmd) && std::filesystem::exists(tmp_path)
                  && std::filesystem::file_size(tmp_path) > 0;
        if (!ok) {
          std::filesystem::remove(tmp_path, ec);
          return false;
        }
        std::filesystem::rename(tmp_path, dest_path, ec);
        if (ec) {
          // `rename` can fail across filesystems; fall back to copy + remove.
          ec.clear();
          std::filesystem::copy_file(tmp_path, dest_path,
                                     std::filesystem::copy_options::overwrite_existing, ec);
          std::filesystem::remove(tmp_path, ec);
        }
        return !ec;
      }

      /// Fetch the contents of `url` (e.g. an HTML directory listing) as a
      /// string. Returns std::nullopt if the request fails or downloads are
      /// disabled.
      std::optional<std::string> FetchTextFromUrl(const std::string& url) {
        if (IsKernelDownloadDisabled()) return std::nullopt;

        auto tmp_path
            = std::filesystem::temp_directory_path()
              / fmt::format("lupnt_naif_listing_{}.html", static_cast<long>(std::time(nullptr)));
        std::string cmd
            = fmt::format("curl -fsSL --connect-timeout 10 --max-time 60 -o \"{}\" \"{}\"",
                          tmp_path.string(), url);
        std::error_code ec;
        if (!RunShellCommand(cmd) || !std::filesystem::exists(tmp_path)) {
          std::filesystem::remove(tmp_path, ec);
          return std::nullopt;
        }
        std::ifstream ifs(tmp_path);
        std::stringstream ss;
        ss << ifs.rdbuf();
        ifs.close();
        std::filesystem::remove(tmp_path, ec);
        return ss.str();
      }

      /// Find the "latest" file name matching `pattern` in the HTML
      /// directory listing at `listing_url`. NAIF kernel file names embed
      /// zero-padded version numbers and dates (e.g.
      /// `moon_pa_de440_200625.bpc`, `moon_de440_250416.tf`), so the
      /// lexicographically-largest match also corresponds to the most
      /// recently produced kernel. Returns std::nullopt if the listing could
      /// not be retrieved or no file name matched.
      std::optional<std::string> FindLatestFilenameInListing(const std::string& listing_url,
                                                             const std::regex& pattern) {
        auto contents = FetchTextFromUrl(listing_url);
        if (!contents.has_value()) return std::nullopt;

        std::string latest;
        for (auto it = std::sregex_iterator(contents->begin(), contents->end(), pattern);
             it != std::sregex_iterator(); ++it) {
          const std::string match = it->str();
          if (match > latest) latest = match;
        }
        if (latest.empty()) return std::nullopt;
        return latest;
      }

      /// Find the "latest" (lexicographically-largest) file name matching
      /// `pattern` already cached in `dir`. Used as an offline fallback when
      /// the NAIF directory listing cannot be retrieved.
      std::optional<std::string> FindLatestFilenameLocally(const std::filesystem::path& dir,
                                                           const std::regex& pattern) {
        std::string latest;
        std::error_code ec;
        for (const auto& entry : std::filesystem::directory_iterator(dir, ec)) {
          if (!entry.is_regular_file()) continue;
          std::string name = entry.path().filename().string();
          if (std::regex_match(name, pattern) && name > latest) latest = name;
        }
        if (latest.empty()) return std::nullopt;
        return latest;
      }

      /// Ensure that the kernel file `filename` is present in `kernel_dir` —
      /// (re)downloading the current copy from `url` when network access is
      /// available — and load it into the kernel pool with `furnsh_c`. If the
      /// download fails (e.g. no network access), the existing local copy (if
      /// any) is loaded instead so that LuPNT keeps working offline. Assumes
      /// the current working directory is `kernel_dir`. Returns true if the
      /// kernel ended up loaded.
      bool RefreshAndLoadKernel(const std::string& filename, const std::string& url,
                                const std::filesystem::path& kernel_dir) {
        std::filesystem::path local_path = kernel_dir / filename;
        bool had_local_copy = std::filesystem::exists(local_path);

        if (DownloadToFile(url, local_path)) {
          Logger::Info(fmt::format("SPICE: downloaded latest NAIF kernel '{}'", filename), "Spice");
        } else if (had_local_copy) {
          Logger::Debug(
              fmt::format("SPICE: could not refresh '{}' from NAIF; using cached local copy",
                          filename),
              "Spice");
        } else {
          Logger::Warn(fmt::format("SPICE: could not download '{}' from NAIF ({}) and no local "
                                   "copy is available; this kernel will not be loaded",
                                   filename, url),
                       "Spice");
          return false;
        }

        furnsh_c(filename.c_str());
        CheckSpiceFailure(fmt::format("Loading {}", filename));
        return true;
      }
    }  // namespace

    /**
     * @brief load the Spice kernels
     *
     */
    void LoadSpiceKernel(void) {
      if (spice_loaded) return;
      SpiceInt kcount;
      ktotal_c("ALL", &kcount);
      if (kcount > 0) return;

      SetSpiceReturnMode();

      std::string orig_dir = std::filesystem::current_path().string();
      std::filesystem::path kernel_dir = GetCspiceKernelDir();
      std::filesystem::current_path(kernel_dir);

      furnsh_c("naif0012.tls");  // leap seconds
      CheckSpiceFailure("Loading naif0012.tls");
      furnsh_c("de440.bsp");  // planetary ephemeris
      CheckSpiceFailure("Loading de440.bsp");

      // High-fidelity TT-TDB time ephemeris. de440t.bsp is de440.bsp plus a
      // Chebyshev segment (NAIF body 1000000001) giving the integrated
      // TT - TDB difference at the geocenter. ConvertTime uses it for TDB<->TT
      // in place of the truncated analytic series in unitim_c (the two differ
      // by up to ~30 us). The extra planetary segments duplicate de440.bsp and
      // are harmless (SPICE gives precedence to the last-loaded file).
      if (std::filesystem::exists("de440t.bsp")) {
        furnsh_c("de440t.bsp");
        CheckSpiceFailure("Loading de440t.bsp");
        tt_tdb_kernel_loaded = true;
      } else {
        Logger::Warn(
            "SPICE: de440t.bsp not found in the kernel directory; ConvertTime "
            "TDB<->TT will fall back to the lower-fidelity analytic series "
            "(unitim_c) instead of the DE440t TT-TDB ephemeris",
            "Spice");
      }

      furnsh_c("pck00011.tpc");  // planetary constants
      CheckSpiceFailure("Loading pck00011.tpc");

      // ----------------------------------------------------------------
      // High-accuracy Earth and lunar orientation kernels
      //
      // The default text PCK ("pck00011.tpc") only provides the low-accuracy
      // IAU_EARTH / IAU_MOON body-fixed frames (errors of roughly 0.06 deg /
      // 6 km for Earth, and up to ~0.005 deg / ~155 m for the Moon).
      // High-accuracy orientation requires binary PCKs that define the
      // ITRF93 (Earth) and MOON_PA_DExxx / MOON_ME_DExxx (Moon) frames; see
      // the NAIF tutorial "'High Accuracy' Orientation and Body-fixed Frames
      // for the Moon and Earth" (April 2023).
      //
      // NAIF continually republishes these kernels as new orientation data
      // becomes available, so during initialization LuPNT looks up and
      // downloads the *current latest* files directly from the NAIF generic
      // kernels server, caching them under the local kernel directory and
      // falling back to the cached copy (or older bundled kernels) if the
      // server cannot be reached. Set LUPNT_SKIP_SPICE_KERNEL_DOWNLOAD=1 to
      // force fully-offline operation using only locally-cached kernels.
      // ----------------------------------------------------------------
      const std::string kPckUrl = std::string(kNaifGenericKernelsUrl) + "/pck/";
      const std::string kLunarFkUrl = std::string(kNaifGenericKernelsUrl) + "/fk/satellites/";

      // --- Earth --------------------------------------------------------
      // NAIF publishes the latest reconstructed + short-term-predicted
      // high-accuracy Earth orientation under a stable file name that never
      // changes ("earth_latest_high_prec.bpc"), making it trivial to always
      // fetch the most recent data.
      const std::string kEarthLatestPck = "earth_latest_high_prec.bpc";
      bool earth_loaded
          = RefreshAndLoadKernel(kEarthLatestPck, kPckUrl + kEarthLatestPck, kernel_dir);
      if (!earth_loaded) {
        // Fall back to whichever (older, hard-coded) Earth orientation
        // kernels may already be cached locally from a previous LuPNT
        // release, so that LuPNT keeps working fully offline.
        for (const char* fallback :
             {"earth_200101_990825_predict.bpc", "earth_000101_241014_240722.bpc"}) {
          if (std::filesystem::exists(fallback)) {
            furnsh_c(fallback);
            CheckSpiceFailure(fmt::format("Loading {}", fallback));
          }
        }
      }

      // --- Moon ---------------------------------------------------------
      // Unlike Earth, NAIF does not publish the latest lunar orientation
      // kernels under a stable file name: a new binary PCK + frame kernel
      // (FK) pair is released -- named after the underlying JPL DE/LE
      // ephemeris model and release date -- each time a new ephemeris is
      // produced (e.g. "moon_pa_de440_200625.bpc" + "moon_de440_250416.tf").
      // Determine the current latest pair by inspecting the NAIF directory
      // listings: file names embed zero-padded DE version numbers and
      // dates, so the lexicographically-largest match is also the most
      // recent. If the listings cannot be retrieved, fall back to the
      // newest matching kernel already cached locally.
      static const std::regex kMoonPckPattern(R"(moon_pa_de\d+_[0-9]{6}\.bpc)");
      static const std::regex kMoonFkPattern(R"(moon_de\d+_[0-9]{6}\.tf)");

      auto moon_pck_name = FindLatestFilenameInListing(kPckUrl, kMoonPckPattern);
      if (!moon_pck_name.has_value())
        moon_pck_name = FindLatestFilenameLocally(kernel_dir, kMoonPckPattern);

      auto moon_fk_name = FindLatestFilenameInListing(kLunarFkUrl, kMoonFkPattern);
      if (!moon_fk_name.has_value())
        moon_fk_name = FindLatestFilenameLocally(kernel_dir, kMoonFkPattern);

      // Load the binary PCK before the frame kernel: the FK connects the
      // MOON_PA_DExxx frame name to the orientation data provided by the PCK.
      if (moon_pck_name.has_value()) {
        RefreshAndLoadKernel(*moon_pck_name, kPckUrl + *moon_pck_name, kernel_dir);
      } else {
        Logger::Warn(
            "SPICE: could not determine the latest lunar orientation binary PCK from "
            "NAIF; high-accuracy lunar orientation (MOON_PA / MOON_ME) will be "
            "unavailable",
            "Spice");
      }
      if (moon_fk_name.has_value()) {
        RefreshAndLoadKernel(*moon_fk_name, kLunarFkUrl + *moon_fk_name, kernel_dir);
      } else {
        Logger::Warn(
            "SPICE: could not determine the latest lunar frame kernel (FK) from NAIF; "
            "high-accuracy lunar orientation (MOON_PA / MOON_ME) will be unavailable",
            "Spice");
      }

      // Make the high-accuracy MOON_PA frame the default lunar body-fixed
      // frame for the legacy (pre-N0062) SPICE APIs (LSPCN, ET2LST, ILLUM,
      // SRFXPT, SUBPT, SUBSOL, ...).
      const std::string kMoonAssocPa = "moon_assoc_pa.tf";
      RefreshAndLoadKernel(kMoonAssocPa, kLunarFkUrl + kMoonAssocPa, kernel_dir);

      // Mars. NAIF ships this as "mar097.bsp"; the guard previously tested for
      // "mars097.bsp", which never matched, so the kernel was silently never
      // furnished. Both spellings are accepted so that either bundle works.
      for (const char* mars_bsp : {"mar097.bsp", "mars097.bsp"}) {
        if (std::filesystem::exists(mars_bsp)) {
          furnsh_c(mars_bsp);
          CheckSpiceFailure(fmt::format("Loading {}", mars_bsp));
          break;
        }
      }

      // Load Chebyshev coefficients
      cheby_s = spk_extract("de440.bsp", &cheby_n);
      if (cheby_s == nullptr) {
        throw std::runtime_error(
            "Could not load SPK file - Please Download the SPK file. See "
            "data/ephemeris/readme.md for instructions");
      }

      std::filesystem::current_path(orig_dir);
      spice_loaded = true;
      return;
    }

    void LoadSpiceKernel(const std::filesystem::path& filepath) {
      LoadSpiceKernel();
      SetSpiceReturnMode();
      furnsh_c(filepath.string().c_str());
      CheckSpiceFailure(fmt::format("Loading {}", filepath.string()));
    }

    int GetNaifId(const std::string& body_name) {
      LoadSpiceKernel();
      SpiceInt code;
      SpiceBoolean found;
#pragma omp critical
      {
        bodn2c_c(body_name.c_str(), &code, &found);
      }
      CheckSpiceFailure(fmt::format("Resolving NAIF body {}", body_name));
      if (!found) throw std::runtime_error(fmt::format("NAIF body not found: {}", body_name));
      return static_cast<int>(code);
    }

    std::string GetNaifName(int naif_id) {
      LoadSpiceKernel();
      SpiceChar name[256];
      SpiceBoolean found;
#pragma omp critical
      {
        bodc2n_c(static_cast<SpiceInt>(naif_id), sizeof(name), name, &found);
      }
      CheckSpiceFailure(fmt::format("Resolving NAIF body id {}", naif_id));
      if (!found) return std::to_string(naif_id);
      return std::string(name);
    }

    bool HasNaifBody(const std::string& body_name) {
      LoadSpiceKernel();
      SpiceInt code;
      SpiceBoolean found;
#pragma omp critical
      {
        bodn2c_c(body_name.c_str(), &code, &found);
      }
      CheckSpiceFailure(fmt::format("Resolving NAIF body {}", body_name));
      return static_cast<bool>(found);
    }

    bool HasFrame(const std::string& frame_name) {
      LoadSpiceKernel();
      SpiceInt frame_code = 0;
#pragma omp critical
      {
        namfrm_c(frame_name.c_str(), &frame_code);
      }
      CheckSpiceFailure(fmt::format("Resolving NAIF frame {}", frame_name));
      return frame_code != 0;
    }

    /**
     * @brief Extracts the PCK coefficients from the kernel files
     *
     */
    void ExtractPckCoeffs() {
      int handle;
      SpiceInt pck_handle;
      ConstSpiceChar* pck_file = "../data/ephemeris/moon_pa_de440_200625.bpc";
      double t_tdb = 8000.0;
      int body = 301;
      double descr[5];
      char ident[40];
      int found;
      SpiceDouble record[120];
      int rsize;
      int pdeg;
      // double ra;
      // double dec;
      // double w;
      // double lambda;
      static char bref[32];
      double eulang[6];

      // The fundamental quantities defined by PCK orientation models are actually
      // Euler angles, not matrices. These Euler angles, which we call ``RA, DEC,
      // and W,'' are related to the transformation operator returned from
      // pxform_c by the equation rotate = [ W ]   [ Pi/2 - DEC ]   [ Pi/2 + RA ]
      //              3                1               3
      // To directly retrieve these angles, use the call:

      // bodeul_( &body, &t_tdb, &ra, &dec, &w, &lambda );
      // std::cout << "body: " << body << std::endl;
      // std::cout << "t_tdb: " << t_tdb << std::endl;
      // std::cout << "ra: " << ra << std::endl;
      // std::cout << "dec: " << dec << std::endl;
      // std::cout << "w: " << w << std::endl;
      // std::cout << " " << std::endl;

      pckeul_(&body, &t_tdb, &found, bref, eulang, (ftnlen)32);
      std::cout << "found:" << found << std::endl;
      std::cout << "phi: " << eulang[0] << std::endl;
      std::cout << "delta: " << eulang[1] << std::endl;
      std::cout << " " << std::endl;

      pcklof_c(pck_file, &pck_handle);  // load the PCK file
      pcksfs_(&body, &t_tdb, &handle, descr, ident, &found, (ftnlen)40);

      std::cout << "pck handle: :" << pck_handle << std::endl;
      std::cout << "handle: :" << handle << std::endl;
      std::cout << "descr: " << &descr << std::endl;
      std::cout << "ident: " << &ident << std::endl;
      std::cout << "found:" << found << std::endl;

      if (found) {
        pckr02_(&handle, descr, &t_tdb, record);
        rsize = record[1];
        pdeg = (rsize - 2) / 3 - 1;
        std::cout << "Polynomial Size:" << rsize << std::endl;
        std::cout << "Polynomial Degree:" << pdeg << std::endl;
      }

      // extract coefficients from the CSPICE PCK file
      //    pcksfs_(body, t_tdb, &handle, descr, ident, found, (ftnlen)40);
      //    pckr02_c(handle, target)
    }

    Vec3d GetBodyPosSpice(Real t_tdb, BodyId obs, BodyId target, const std::string& refFrame,
                          const std::string& abCorrection) {
      if (!spice_loaded) {
        LoadSpiceKernel();
      }

      std::string targ_str = std::to_string((int)target);
      std::string obs_str = std::to_string((int)obs);

      SpiceDouble ptarg[3];
      SpiceDouble et = t_tdb.val();
      const char* targ = strcpy(new char[targ_str.length() + 1], targ_str.c_str());
      const char* ref = strcpy(new char[refFrame.length() + 1], refFrame.c_str());
      const char* abcorr = strcpy(new char[abCorrection.length() + 1], abCorrection.c_str());
      const char* obs_spice = strcpy(new char[obs_str.length() + 1], obs_str.c_str());
      SpiceDouble lt;

      // void spkpos_c(ConstSpiceChar * targ, SpiceDouble t_tdb, ConstSpiceChar *
      // ref, ConstSpiceChar * abcorr,
      //               ConstSpiceChar * obs, SpiceDouble ptarg[3], SpiceDouble *
      //               lt)
      spkpos_c(targ, et, ref, abcorr, obs_spice, ptarg, &lt);

      Vec3d r;
      for (int i = 0; i < 3; i++) r(i) = ptarg[i];
      return r;
    }

    Vec6d GetBodyPosVelSpice(Real t_tdb, BodyId obs, BodyId target, const std::string& refFrame,
                             const std::string& abCorrection) {
      if (!spice_loaded) LoadSpiceKernel();

      SpiceDouble starg[6];
      Vec6d rv;

      std::string targ_str = std::to_string((int)target);
      std::string obs_str = std::to_string((int)obs);

      SpiceDouble et = t_tdb.val();
      const char* targ = strcpy(new char[targ_str.length() + 1], targ_str.c_str());
      const char* ref = strcpy(new char[refFrame.length() + 1], refFrame.c_str());
      const char* abcorr = strcpy(new char[abCorrection.length() + 1], abCorrection.c_str());
      const char* obs_spice = strcpy(new char[obs_str.length() + 1], obs_str.c_str());
      SpiceDouble lt;

      //  void spkez_c ( SpiceInt            targ,
      //                 SpiceDouble         et,
      //                 ConstSpiceChar     *ref,
      //                 ConstSpiceChar     *abcorr,
      //                 SpiceInt            obs,
      //                 SpiceDouble         starg[6],
      //                 SpiceDouble        *lt        )
#pragma omp critical
      {
        spkezr_c(targ, et, ref, abcorr, obs_spice, starg, &lt);
      }
      for (int i = 0; i < 6; i++) rv(i) = starg[i];
      return rv;
    }

    GroundStationSpiceData GetGroundStationDataSpice(Real t_tdb, const std::string& station_name,
                                                     BodyId center, const std::string& refFrame,
                                                     const std::string& abCorrection) {
      int station_id = GetNaifId(station_name);
      GroundStationSpiceData data
          = GetGroundStationDataSpice(t_tdb, station_id, center, refFrame, abCorrection);
      data.name = station_name;
      return data;
    }

    GroundStationSpiceData GetGroundStationDataSpice(Real t_tdb, int station_id, BodyId center,
                                                     const std::string& refFrame,
                                                     const std::string& abCorrection) {
      LoadSpiceKernel();
      SetSpiceReturnMode();

      SpiceDouble et = t_tdb.val();
      SpiceDouble pos_km[3];
      SpiceDouble lt;
      std::string station = std::to_string(station_id);
      std::string center_id = std::to_string(static_cast<int>(center));

#pragma omp critical
      {
        spkpos_c(station.c_str(), et, refFrame.c_str(), abCorrection.c_str(), center_id.c_str(),
                 pos_km, &lt);
      }
      CheckSpiceFailure(fmt::format("Fetching ground station {} from SPICE", station_id));

      GroundStationSpiceData data;
      data.naif_id = station_id;
      data.name = GetNaifName(station_id);
      data.body_id = center;
      data.frame = refFrame;
      for (int i = 0; i < 3; ++i) data.position_m(i) = pos_km[i] * M_KM;

      BodyData body_data = GetBodyData(center);
      Cart3 pos(data.position_m.cast<Real>(), body_data.fixed_frame);
      LatLonAlt lla = CartToLatLonAlt(pos, body_data.R, body_data.flattening);
      data.latitude_deg = lla(0).val();
      data.longitude_deg = lla(1).val();
      data.altitude_m = lla(2).val();
      return data;
    }

    /**
     * @brief Get the Frame Conversion Mat object
     *
     * @param t_tdb
     * @param from_frame
     * @param to_frame
     * @return VecXd
     */
    Mat6d GetFrameConversionMat(Real t_tdb, const std::string& from_frame,
                                const std::string& to_frame) {
      if (!spice_loaded) LoadSpiceKernel();

      SpiceDouble et_spice = (SpiceDouble)t_tdb.val();
      double xform[6][6];
      Mat6d M_rot;

      const char* from_frame_char = strcpy(new char[from_frame.length() + 1], from_frame.c_str());
      const char* to_frame_char = strcpy(new char[to_frame.length() + 1], to_frame.c_str());
#pragma omp critical
      {
        sxform_c(from_frame_char, to_frame_char, et_spice, xform);
      }

      for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
          M_rot(i, j) = xform[i][j];
        }
      }

      return M_rot;
    }

    Vec3d GetPlanetOrientation(BodyId id, Real t_tdb) {
      if (!spice_loaded) LoadSpiceKernel();

      std::string to_frame = "J2000";

      char bodyname[36];
      SpiceBoolean found;
#pragma omp critical
      {
        bodc2n_c((SpiceInt)id, 36, bodyname, &found);
      }
      if (!found) {
        throw std::runtime_error("Invalid planet ID");
      }
      std::string from_frame = "IAU_" + std::string(bodyname);

      SpiceDouble rotmat[3][3];
      SpiceDouble et_spice = (SpiceDouble)t_tdb.val();
      const char* from_frame_char = strcpy(new char[from_frame.length() + 1], from_frame.c_str());
      const char* to_frame_char = strcpy(new char[to_frame.length() + 1], to_frame.c_str());
#pragma omp critical
      {
        pxform_c(from_frame_char, to_frame_char, et_spice, rotmat);
      }
      SpiceDouble psi, theta, phi;

#pragma omp critical
      m2eul_c(rotmat, 3, 1, 3, &psi, &theta, &phi);

      Vec3d angles;
      double W = -(phi + PI);          // R_z(-W)
      double delta0 = PI / 2 - theta;  // R_x(PI/2 - delta0)
      double alpha0 = PI / 2 - psi;    // R_z(-alpha0)
      angles << alpha0, delta0, W;
      return angles;
    }

    /**
     * @brief Convert a string to ephemeris time
     *
     * @param str           string to be converted
     * SO (T) Formats.

       String                        Year Mon  DOY DOM  HR Min Sec
       ----------------------------  ---- ---  --- ---  -- --- ------
       1996-12-18T12:28:28           1996 Dec   na  18  12  28 28
       1986-01-18T12                 1986 Jan   na  18  12  00 00
       1986-01-18T12:19              1986 Jan   na  18  12  19 00
       1986-01-18T12:19:52.18        1986 Jan   na  18  12  19 52.18
       1986-01-18T12:19:52.18Z       1986 Jan   na  18  12  19 52.18
       1995-08T18:28:12              1995  na  008  na  18  28 12
       1995-08T18:28:12Z             1995  na  008  na  18  28 12
       1995-18T                      1995  na  018  na  00  00 00
       0000-01-01T                   1 BC Jan   na  01  00  00 00
    Calendar Formats.
       String                        Year   Mon DOM  HR Min  Sec
       ----------------------------  ----   --- ---  -- ---  ------
       Tue Aug  6 11:10:57  1996     1996   Aug  06  11  10  57
       1 DEC 1997 12:28:29.192       1997   Dec  01  12  28  29.192
       2/3/1996 17:18:12.002         1996   Feb  03  17  18  12.002
       Mar 2 12:18:17.287 1993       1993   Mar  02  12  18  17.287
       1992 11:18:28  3 Jul          1992   Jul  03  11  18  28
       June 12, 1989 01:21           1989   Jun  12  01  21  00
       1978/3/12 23:28:59.29         1978   Mar  12  23  28  59.29
       17JUN1982 18:28:28            1982   Jun  17  18  28  28
       13:28:28.128 1992 27 Jun      1992   Jun  27  13  28  28.128
       1972 27 jun 12:29             1972   Jun  27  12  29  00
       '93 Jan 23 12:29:47.289       1993*  Jan  23  12  29  47.289
       27 Jan 3, 19:12:28.182        2027*  Jan  03  19  12  28.182
       23 A.D. APR 4, 18:28:29.29    0023** Apr  04  18  28  29.29
       18 B.C. Jun 3, 12:29:28.291   -017** Jun  03  12  29  28.291
       29 Jun  30 12:29:29.298       2029+  Jun  30  12  29  29.298
       29 Jun '30 12:29:29.298       2030*  Jun  29  12  29  29.298
    Day of Year Formats.
       String                        Year  DOY HR Min Sec
       ----------------------------  ----  --- -- --- ------
       1997-162::12:18:28.827        1997  162 12  18 28.827
       162-1996/12:28:28.287         1996  162 12  28 28.287
       1993-321/12:28:28.287         1993  231 12  28 28.287
       1992 183// 12:18:19           1992  183 12  18 19
       17:28:01.287 1992-272//       1992  272 17  28 01.287
       17:28:01.282 272-1994//       1994  272 17  28 01.282
       '92-271/ 12:28:30.291         1992* 271 12  28 30.291
       92-182/ 18:28:28.281          1992* 182 18  28 28.281
       182-92/ 12:29:29.192          0182+ 092 12  29 29.192
       182-'92/ 12:28:29.182         1992  182 12  28 29.182
    Julian Date Strings.
       jd 28272.291                  Julian Date   28272.291
       2451515.2981 (JD)             Julian Date 2451515.2981
       2451515.2981 JD               Julian Date 2451515.2981

     * @return real     ephemeris time (TDB) (seconds past the J2000 epoch)
     */
    Real StringToTdb(const std::string& str) {
      if (!spice_loaded) LoadSpiceKernel();
      SpiceDouble t_tdb;
      str2et_c(str.c_str(), &t_tdb);
      return t_tdb;
    }

    /**
     * @brief Convert string to TAI
     *
     * @param str time string
     * @return real
     */
    Real StringToTai(const std::string& str) {
      if (!spice_loaded) LoadSpiceKernel();
      Real t_tdb = StringToTdb(str);
      Real t_tai = spice::ConvertTime(t_tdb, Time::TDB, Time::TAI);
      return t_tai;
    }

    /**
     * @brief Convert string to UTC
     *
     * @param t_tdb time in TDB
     * @param prec precision of the output string (default 3)
     * @return std::string
     */
    std::string TDBtoStringUTC(Real t_tdb, int prec = 3) {
      if (!spice_loaded) LoadSpiceKernel();
      SpiceDouble et = t_tdb.val();
      SpiceChar str[100];
      et2utc_c(et, "C", prec, 100, str);
      return std::string(str);
    }

    /**
     * @brief Convert TAI to string UTC
     *
     * @param t_tai time in TAI
     * @param prec precision of the output string (default 3)
     * @return std::string
     */
    std::string TAItoStringUTC(Real t_tai, int prec = 3) {
      if (!spice_loaded) LoadSpiceKernel();
      Real et_tdb = spice::ConvertTime(t_tai, Time::TAI, Time::TDB);
      std::string str = TDBtoStringUTC(et_tdb, prec);
      return str;
    }

    /**
     * @brief Convert time from one time system to another
     *
     * @param t     in time in seconds
     * @param from  from time system
     *  String ID   Time system
     *  ---------   --------------------------
       TAI         International Atomic Time
       TDB         Barycentric Dynamical Time
       TT          Terrestrial Time
       TDT         Terrestrial Dynamical Time (TT)
       ET          Ephemeris time, alias for TDB
       JDTDB       Julian Date relative to TDB
       JDTDT       Julian Date relative to TDT (TT)
       JED         Julian Ephemeris date (synonym to JDTDB)
       GPS         Global Positioning System Time

     * @param to  to time system
     * @return real     out time in seconds
     */
    /// @brief Read TT - TDB [s] at a TDB epoch from the DE440t time ephemeris.
    ///
    /// The difference is stored as the x-component of the position of NAIF
    /// body 1000000001 ("TT-TDB") relative to 1000000000. The ephemeris is
    /// parameterized by TDB; not thread-safe, so callers guard `spkgeo_c` with
    /// an omp-critical section.
    static double TtMinusTdbFromKernel(double et_tdb) {
      SpiceDouble state[6];
      SpiceDouble lt;
      spkgeo_c(kNaifTtMinusTdb, et_tdb, "J2000", kNaifTdb, state, &lt);
      CheckSpiceFailure("Reading DE440t TT-TDB ephemeris");
      return state[0];  // TT - TDB [s]
    }

    Real GetTimeEphemerisOffset(Real t_tdb, int naif_id) {
      if (!spice_loaded) LoadSpiceKernel();
      double offset = 0.0;
      bool failed = false;
      // Capture the failure flag rather than throwing inside the critical
      // region: an exception must not escape an OpenMP structured block.
#pragma omp critical
      {
        SpiceDouble state[6];
        SpiceDouble lt;
        spkgeo_c(naif_id, t_tdb.val(), "J2000", kNaifTdb, state, &lt);
        failed = failed_c();
        if (failed) {
          reset_c();
        } else {
          offset = state[0];
        }
      }
      LUPNT_CHECK(!failed,
                  fmt::format("Could not read time ephemeris for NAIF id {} at t_tdb={} "
                              "(is the containing kernel loaded?)",
                              naif_id, t_tdb.val()),
                  "Spice");
      return offset;
    }

    Real ConvertTime(Real t, Time from, Time to) {
      if (!spice_loaded) LoadSpiceKernel();
      if (from == to) return t;

      // High-fidelity TT<->TDB via the DE440t TT-TDB ephemeris. unitim_c uses
      // a truncated analytic series (Fairhead & Bretagnon) for this leg; the
      // integrated DE440t difference is more accurate by up to ~30 us. Route
      // any TDB-involving conversion through TT so the TT<->TDB leg uses the
      // kernel while the remaining (leap-second / linear) legs use unitim_c.
      if (tt_tdb_kernel_loaded && (from == Time::TDB || to == Time::TDB) && from != to) {
        if (from == Time::TDB && to == Time::TT) {
          double offset;
#pragma omp critical
          {
            offset = TtMinusTdbFromKernel(t.val());
          }
          return t + offset;  // t is TDB -> TT
        }
        if (from == Time::TT && to == Time::TDB) {
          // The ephemeris is parameterized by TDB; using the TT epoch as the
          // lookup argument adds < 1e-12 s error (|TT-TDB| < 2 ms, rate ~3e-10).
          double offset;
#pragma omp critical
          {
            offset = TtMinusTdbFromKernel(t.val());
          }
          return t - offset;  // t is TT -> TDB
        }
        if (from == Time::TDB)
          return spice::ConvertTime(spice::ConvertTime(t, Time::TDB, Time::TT), Time::TT, to);
        if (to == Time::TDB)
          return spice::ConvertTime(spice::ConvertTime(t, from, Time::TT), Time::TT, Time::TDB);
      }

      SpiceDouble t_in = t.val();
      SpiceDouble t_out_spice;
#pragma omp critical
      {
        t_out_spice
            = unitim_c(t_in, time_to_string.at(from).c_str(), time_to_string.at(to).c_str());
      }
      double offset = t_out_spice - t_in;  // offset in seconds
      Real t_out = t + offset;             // this is to convert to real
      return t_out;
    }

    MatX6 GetBodyPosVel(const VecX& t_tdb, BodyId center, BodyId target) {
      if (!spice_loaded) LoadSpiceKernel();
      MatX6 retState(t_tdb.size(), 6);
      for (int i = 0; i < t_tdb.size(); i++)
        retState.row(i) = GetBodyPosVel(t_tdb(i), center, target).transpose();
      return retState;
    }

    Vec6 GetBodyPosVelBase(const Real t_tdb, BodyId center, BodyId target) {
      if (!spice_loaded) LoadSpiceKernel();
      if (center == target) return Vec6::Zero();

      for (int i = 0; i < cheby_n; i++) {
        if (cheby_s[i].center == (int)center && cheby_s[i].target == (int)target)
          return cheby_posvel_ad(t_tdb, cheby_s[i].seg, cheby_s[i].len);
        if (cheby_s[i].center == (int)target && cheby_s[i].target == (int)center)
          return -cheby_posvel_ad(t_tdb, cheby_s[i].seg, cheby_s[i].len);
      }
      LUPNT_CHECK(false, "Chebyshev coefficients not found", "SpiceInterface");
      return Vec6::Zero();
    }

    /**
     * @brief Get the Body Position and Velocity using Chebyshev polynomials
     *
     * @param t_tdb TDB seconds from J2000.
     * @param center  center body id
     * @param target  target body id
     * @return Vec6 in intertial axes
     */
    Vec6 GetBodyPosVel(const Real t_tdb, BodyId center, BodyId target) {
      if (!spice_loaded) LoadSpiceKernel();
      if (center == target) return Vec6::Zero();

      Vec6 rv = Vec6::Zero();
      if (target == BodyId::EARTH) {
        rv += GetBodyPosVelBase(t_tdb, BodyId::EMB, BodyId::EARTH);
        target = BodyId::EMB;
      }
      if (center == BodyId::EARTH) {
        rv -= GetBodyPosVelBase(t_tdb, BodyId::EMB, BodyId::EARTH);
        center = BodyId::EMB;
      }
      if (target == BodyId::MOON) {
        rv += GetBodyPosVelBase(t_tdb, BodyId::EMB, BodyId::MOON);
        target = BodyId::EMB;
      }
      if (center == BodyId::MOON) {
        rv -= GetBodyPosVelBase(t_tdb, BodyId::EMB, BodyId::MOON);
        center = BodyId::EMB;
      }
      if (target == BodyId::MERCURY) {
        rv += GetBodyPosVelBase(t_tdb, BodyId::MERCURY_BARYCENTER, BodyId::MERCURY);
        target = BodyId::MERCURY_BARYCENTER;
      }
      if (center == BodyId::MERCURY) {
        rv -= GetBodyPosVelBase(t_tdb, BodyId::MERCURY_BARYCENTER, BodyId::MERCURY);
        center = BodyId::MERCURY_BARYCENTER;
      }
      if (target == BodyId::VENUS) {
        rv += GetBodyPosVelBase(t_tdb, BodyId::VENUS_BARYCENTER, BodyId::VENUS);
        target = BodyId::VENUS_BARYCENTER;
      }
      if (center == BodyId::VENUS) {
        rv -= GetBodyPosVelBase(t_tdb, BodyId::VENUS_BARYCENTER, BodyId::VENUS);
        center = BodyId::VENUS_BARYCENTER;
      }
      rv += GetBodyPosVelBase(t_tdb, center, target);
      return rv;
    }

  }  // namespace spice

}  // namespace lupnt
