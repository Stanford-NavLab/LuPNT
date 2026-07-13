#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief One high-resolution lunar-south-pole landing site published in NASA PGDA
  /// product 78 ("High-Resolution LOLA Topography for Lunar South Pole Sites").
  ///
  /// The `lat_deg`/`lon_deg` are the *approximate* geographic center of the site tile,
  /// used only to pick the nearest site for a query latitude/longitude (see
  /// `SelectLolaSite`). The authoritative extent/coordinates come from the GeoTIFF's own
  /// geotransform when the tile is loaded.
  struct LolaSite {
    std::string id;    ///< PGDA folder id, e.g. "Site01".
    std::string name;  ///< Human-readable name, e.g. "Connecting ridge".
    double lat_deg;    ///< Approximate center latitude [deg] (south pole ~ -90).
    double lon_deg;    ///< Approximate center east longitude [deg], 0-360.
  };

  /// @brief The six PGDA product-78 LOLA 5 m/pixel south-pole DEM sites currently hosted
  /// under `https://pgda.gsfc.nasa.gov/data/LOLA_5mpp/`.
  const std::vector<LolaSite>& GetLolaSites();

  /// @brief Pick the PGDA product-78 site whose center is angularly closest to the query
  /// latitude/longitude (great-circle distance on the unit sphere, so it is well behaved
  /// right at the pole where longitude is ill-defined).
  ///
  /// @param lat_deg Query latitude [deg] (negative near the south pole).
  /// @param lon_deg Query east longitude [deg].
  /// @return        The nearest `LolaSite` from `GetLolaSites()`.
  const LolaSite& SelectLolaSite(double lat_deg, double lon_deg);

  /// @brief Build the PGDA product-78 download URL for a site's 5 m/pixel surface DEM
  /// GeoTIFF (`<site>/<site>_final_adj_5mpp_surf.tif`).
  std::string LolaDemUrl(const std::string& site_id);

  /// @brief Download (and cache) a site's 5 m/pixel surface DEM GeoTIFF from NASA PGDA.
  ///
  /// The file is cached at `GetDataPath()/dem/LOLA_5mpp/<site>/<site>_final_adj_5mpp_surf.tif`;
  /// if that already exists it is returned without any network access. Downloads use the
  /// same `curl -fsSL` idiom as the other lupnt loaders (`sp3_loader`, `eop`, `spice`) and
  /// are disabled when the environment variable `LUPNT_SKIP_DEM_DOWNLOAD` is set truthy
  /// (in which case a missing cache throws).
  ///
  /// @param site_id PGDA site folder id (e.g. "Site01").
  /// @return        Local filesystem path to the cached GeoTIFF.
  std::filesystem::path DownloadLolaDem(const std::string& site_id);

  /// @brief A loaded (cropped, downsampled) lunar digital elevation model.
  ///
  /// Holds a regular grid of terrain elevations together with the pixel-center x/y
  /// coordinate grids in the GeoTIFF's *native* projected meters (south-polar
  /// stereographic, MOON_ME frame for PGDA product 78). Elevation lookups are bilinear in
  /// those native coordinates; a latitude/longitude convenience lookup is also provided
  /// (approximate -- see `GetElevationLatLon`).
  class LunarDem {
  public:
    LunarDem() = default;

    /// @brief Construct from LoadTiff-style outputs plus the originating site.
    /// @param x    Pixel-center x-coordinate grid [m] (native projection), shape (ny, nx).
    /// @param y    Pixel-center y-coordinate grid [m] (native projection), shape (ny, nx).
    /// @param elev Elevation grid [m], shape (ny, nx).
    /// @param site Originating PGDA site metadata.
    LunarDem(const MatXd& x, const MatXd& y, const MatXd& elev, const LolaSite& site);

    /// @brief Bilinearly interpolate the terrain elevation at native projected (x, y).
    /// Queries outside the grid are clamped to the nearest edge.
    /// @param x Native x-coordinate [m].
    /// @param y Native y-coordinate [m].
    /// @return  Elevation [m].
    double GetElevation(double x, double y) const;

    /// @brief Elevation at a geographic latitude/longitude, mapping through lupnt's
    /// polar-stereographic conversion (`LatLonAltToStereographic`, R = R_MOON) before the
    /// bilinear lookup. Approximate: lupnt's generic conformal stereographic is not
    /// guaranteed to match the GeoTIFF's exact projection scale, so prefer `GetElevation`
    /// in native coordinates on the critical path.
    /// @param lat_deg Latitude [deg]. @param lon_deg East longitude [deg]. @return Elevation [m].
    double GetElevationLatLon(double lat_deg, double lon_deg) const;

    const MatXd& x() const { return x_; }             ///< x-coordinate grid [m].
    const MatXd& y() const { return y_; }             ///< y-coordinate grid [m].
    const MatXd& elevation() const { return elev_; }  ///< elevation grid [m].
    const LolaSite& site() const { return site_; }    ///< originating site.

    int rows() const { return static_cast<int>(elev_.rows()); }  ///< grid height (ny).
    int cols() const { return static_cast<int>(elev_.cols()); }  ///< grid width (nx).

    double x_min() const { return std::min(x_left_, x_right_); }    ///< min native x [m].
    double x_max() const { return std::max(x_left_, x_right_); }    ///< max native x [m].
    double y_min() const { return std::min(y_top_, y_bottom_); }    ///< min native y [m].
    double y_max() const { return std::max(y_top_, y_bottom_); }    ///< max native y [m].
    double center_x() const { return 0.5 * (x_left_ + x_right_); }  ///< native x center [m].
    double center_y() const { return 0.5 * (y_top_ + y_bottom_); }  ///< native y center [m].

  private:
    MatXd x_, y_, elev_;
    LolaSite site_;

    // Regular-grid parameters cached for fast bilinear lookup. Columns increase in x,
    // rows increase in y (dy_ is typically negative for a north-up raster).
    double x0_ = 0.0, y0_ = 0.0;  // coordinate of pixel-center (row 0, col 0)
    double dx_ = 1.0, dy_ = 1.0;  // per-column / per-row coordinate step [m]
    int nx_ = 0, ny_ = 0;
    double x_left_ = 0.0, x_right_ = 0.0, y_top_ = 0.0, y_bottom_ = 0.0;
  };

  /// @brief Load the appropriate PGDA product-78 DEM for a query latitude/longitude.
  ///
  /// Selects the nearest site (`SelectLolaSite`), ensures its GeoTIFF is downloaded/cached
  /// (`DownloadLolaDem`), then crops a `2*half_width_m` square window centered on the tile
  /// and downsamples it to at most `max_res` grid spacing via the existing `LoadTiff`.
  ///
  /// When `dem_file` is non-empty it overrides the site download entirely: the given GeoTIFF
  /// is loaded/cropped/downsampled directly (no `DownloadLolaDem`, no network access). The
  /// site metadata attached to the returned `LunarDem` is still the nearest `SelectLolaSite`
  /// entry for the query lat/lon (a synthetic label, since the explicit file may be any tile).
  /// This is the path used to run surface scenarios against a small bundled fixture tile.
  ///
  /// @param lat_deg      Query latitude [deg] (negative near the south pole).
  /// @param lon_deg      Query east longitude [deg].
  /// @param half_width_m Half-width of the square crop window [m] about the tile center.
  /// @param max_res      Maximum output grid spacing [m] (downsampling target).
  /// @param dem_file     Optional explicit GeoTIFF path; if set, skips site select/download.
  /// @return             The cropped/downsampled `LunarDem`.
  LunarDem LoadLolaDem(double lat_deg, double lon_deg, double half_width_m = 5000.0,
                       double max_res = 20.0, const std::filesystem::path& dem_file = {});

}  // namespace lupnt
