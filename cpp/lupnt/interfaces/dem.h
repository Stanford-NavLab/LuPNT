#pragma once

#include <cpl_conv.h>
#include <gdal_priv.h>
#include <gdalwarper.h>

#include <array>
#include <filesystem>
#include <tuple>

#include "lupnt/core/definitions.h"

namespace lupnt {

  /// @brief Load a (possibly cropped and downsampled) digital elevation model from a GeoTIFF file
  /// via GDAL.
  ///
  /// Used by surface/rover environment and dynamics models to obtain lunar terrain elevation
  /// data over a region of interest for terrain-relative navigation and surface dynamics. Opens
  /// `path` with GDAL, intersects the requested x/y bounding box with the raster's extent,
  /// downsamples (nearest-neighbor) so the output grid spacing does not exceed `max_res`, and
  /// returns the cropped/downsampled elevation grid together with its x/y coordinate grids.
  ///
  /// @param path   Path to the GeoTIFF (DEM) file to load.
  /// @param xlims  Desired x-coordinate (e.g. easting) bounds [min, max] in the raster's
  ///               coordinate units (m), intersected with the raster's actual extent.
  /// @param ylims  Desired y-coordinate (e.g. northing) bounds [min, max] in the raster's
  ///               coordinate units (m), intersected with the raster's actual extent.
  /// @param max_res Maximum output grid spacing [m]; the raster is downsampled by an integer
  ///                factor so neither axis' spacing exceeds this value.
  /// @return       Tuple of (x_coords, y_coords, data): pixel-center x and y coordinate grids
  ///               [m] and the corresponding elevation/raster values, all as matrices of shape
  ///               (out_height, out_width).
  std::tuple<MatXd, MatXd, MatXd> LoadTiff(std::filesystem::path path, std::array<double, 2> xlims,
                                           std::array<double, 2> ylims, double max_res);
}  // namespace lupnt
