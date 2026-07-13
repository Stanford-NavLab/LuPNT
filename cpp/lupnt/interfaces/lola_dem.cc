#include "lupnt/interfaces/lola_dem.h"

#include <cpl_conv.h>
#include <fmt/format.h>
#include <gdal_priv.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <system_error>

#include "lupnt/conversions/coordinate_conversions.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/error.h"
#include "lupnt/core/file.h"
#include "lupnt/interfaces/dem.h"

namespace lupnt {

  namespace {
    constexpr char kLolaBaseUrl[] = "https://pgda.gsfc.nasa.gov/data/LOLA_5mpp";

    bool IsTruthyEnv(const char* value) {
      if (value == nullptr) return false;
      std::string s(value);
      std::transform(s.begin(), s.end(), s.begin(),
                     [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      return !(s.empty() || s == "0" || s == "false" || s == "off" || s == "no");
    }

    bool IsDemDownloadDisabled() { return IsTruthyEnv(std::getenv("LUPNT_SKIP_DEM_DOWNLOAD")); }

    /// Download `url` to `dest_path` using curl, writing to a temporary file first and
    /// atomically moving it into place on success (mirrors `data/eop.cc`).
    bool DownloadFile(const std::string& url, const std::filesystem::path& dest_path) {
      std::filesystem::path tmp_path = dest_path;
      tmp_path += ".part";
      std::error_code ec;
      std::filesystem::remove(tmp_path, ec);
      std::filesystem::create_directories(dest_path.parent_path(), ec);

      std::string cmd
          = fmt::format("curl -fsSL --connect-timeout 10 --max-time 1800 -o \"{}\" \"{}\"",
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

    /// Unit direction vector (Moon-fixed) of a geographic latitude/longitude, used for a
    /// pole-safe angular distance between sites.
    std::array<double, 3> UnitVec(double lat_deg, double lon_deg) {
      double lat = lat_deg * RAD, lon = lon_deg * RAD;
      return {std::cos(lat) * std::cos(lon), std::cos(lat) * std::sin(lon), std::sin(lat)};
    }
  }  // namespace

  const std::vector<LolaSite>& GetLolaSites() {
    // Approximate site centers (deg). Only used to select the nearest tile; the exact
    // extent is read from each GeoTIFF's geotransform when loaded. Coordinates follow the
    // PGDA product-78 site list (Connecting ridge, Shackleton, de Gerlache, etc.).
    static const std::vector<LolaSite> kSites = {
        {"Site01", "Connecting ridge", -89.45, 222.8},
        {"Site04", "Shackleton rim", -89.68, 129.2},
        {"Site07", "Peak near Shackleton", -89.42, 137.0},
        {"Site11", "de Gerlache rim", -88.71, 290.0},
        {"Site20", "Leibnitz beta plateau", -85.33, 33.0},
        {"Site23", "Malapert massif", -85.99, 2.9},
    };
    return kSites;
  }

  const LolaSite& SelectLolaSite(double lat_deg, double lon_deg) {
    const auto& sites = GetLolaSites();
    auto q = UnitVec(lat_deg, lon_deg);
    const LolaSite* best = &sites.front();
    double best_dot = -2.0;
    for (const auto& s : sites) {
      auto v = UnitVec(s.lat_deg, s.lon_deg);
      double dot = q[0] * v[0] + q[1] * v[1] + q[2] * v[2];  // larger = closer
      if (dot > best_dot) {
        best_dot = dot;
        best = &s;
      }
    }
    return *best;
  }

  std::string LolaDemUrl(const std::string& site_id) {
    return fmt::format("{}/{}/{}_final_adj_5mpp_surf.tif", kLolaBaseUrl, site_id, site_id);
  }

  std::filesystem::path DownloadLolaDem(const std::string& site_id) {
    std::filesystem::path dest
        = GetDataPath() / "dem" / "LOLA_5mpp" / site_id / (site_id + "_final_adj_5mpp_surf.tif");
    std::error_code ec;
    if (std::filesystem::exists(dest, ec) && std::filesystem::file_size(dest, ec) > 0) return dest;

    LUPNT_CHECK(!IsDemDownloadDisabled(),
                "LOLA DEM not cached and downloads are disabled (LUPNT_SKIP_DEM_DOWNLOAD).",
                "DownloadLolaDem");
    bool ok = DownloadFile(LolaDemUrl(site_id), dest);
    LUPNT_CHECK(ok, "Failed to download LOLA DEM from NASA PGDA (product 78).", "DownloadLolaDem");
    return dest;
  }

  LunarDem::LunarDem(const MatXd& x, const MatXd& y, const MatXd& elev, const LolaSite& site)
      : x_(x), y_(y), elev_(elev), site_(site) {
    ny_ = static_cast<int>(elev_.rows());
    nx_ = static_cast<int>(elev_.cols());
    LUPNT_CHECK(nx_ >= 2 && ny_ >= 2, "DEM grid must be at least 2x2", "LunarDem");

    x0_ = x_(0, 0);
    y0_ = y_(0, 0);
    dx_ = x_(0, 1) - x_(0, 0);
    dy_ = y_(1, 0) - y_(0, 0);
    x_left_ = x_(0, 0);
    x_right_ = x_(0, nx_ - 1);
    y_top_ = y_(0, 0);
    y_bottom_ = y_(ny_ - 1, 0);
  }

  double LunarDem::GetElevation(double x, double y) const {
    LUPNT_CHECK(nx_ >= 2 && ny_ >= 2, "DEM not loaded", "LunarDem::GetElevation");
    // Fractional grid indices, clamped so out-of-range queries hit the nearest edge cell.
    double fc = (x - x0_) / dx_;
    double fr = (y - y0_) / dy_;
    fc = std::clamp(fc, 0.0, static_cast<double>(nx_ - 1));
    fr = std::clamp(fr, 0.0, static_cast<double>(ny_ - 1));

    int c0 = std::min(static_cast<int>(std::floor(fc)), nx_ - 2);
    int r0 = std::min(static_cast<int>(std::floor(fr)), ny_ - 2);
    double tc = fc - c0;
    double tr = fr - r0;

    double e00 = elev_(r0, c0);
    double e01 = elev_(r0, c0 + 1);
    double e10 = elev_(r0 + 1, c0);
    double e11 = elev_(r0 + 1, c0 + 1);
    double top = e00 * (1 - tc) + e01 * tc;
    double bot = e10 * (1 - tc) + e11 * tc;
    return top * (1 - tr) + bot * tr;
  }

  double LunarDem::GetElevationLatLon(double lat_deg, double lon_deg) const {
    Vec3 xya = LatLonAltToStereographic(Vec3(lat_deg, lon_deg, 0.0), R_MOON);
    return GetElevation(xya(0).val(), xya(1).val());
  }

  LunarDem LoadLolaDem(double lat_deg, double lon_deg, double half_width_m, double max_res,
                       const std::filesystem::path& dem_file) {
    const LolaSite& site = SelectLolaSite(lat_deg, lon_deg);
    // An explicit `dem_file` overrides the site download entirely (no network access): the
    // given GeoTIFF is cropped/downsampled directly. Otherwise fall back to the cached/
    // downloaded PGDA tile for the nearest site.
    std::filesystem::path path;
    if (!dem_file.empty()) {
      std::error_code ec;
      LUPNT_CHECK(std::filesystem::exists(dem_file, ec),
                  "Explicit LOLA DEM file not found: " + dem_file.string(), "LoadLolaDem");
      path = dem_file;
    } else {
      path = DownloadLolaDem(site.id);
    }

    // Read the raster geotransform to find the tile center (each site tile is centered on
    // its landing site), then crop a square window around it in native projected meters.
    GDALAllRegister();
    GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(path.string().c_str(), GA_ReadOnly));
    LUPNT_CHECK(ds, "Failed to open cached LOLA DEM GeoTIFF", "LoadLolaDem");
    double gt[6];
    ds->GetGeoTransform(gt);
    int w = ds->GetRasterXSize();
    int h = ds->GetRasterYSize();
    double cx = gt[0] + 0.5 * w * gt[1];
    double cy = gt[3] + 0.5 * h * gt[5];
    GDALClose(ds);

    auto [x, y, data] = LoadTiff(path, {cx - half_width_m, cx + half_width_m},
                                 {cy - half_width_m, cy + half_width_m}, max_res);
    return LunarDem(x, y, data, site);
  }

}  // namespace lupnt
