#include "lupnt/data/dem.h"

#include <cmath>
#include <vector>

#include "lupnt/core/error.h"

namespace lupnt {
  /// @brief Load a TIFF file and return the data, x coordinates, and y coordinates
  /// @param tif_path Path to the TIFF file
  /// @param xlims X limits in m
  /// @param ylims Y limits in m
  /// @param max_res Maximum resolution in m
  /// @return Tuple containing the data, x coordinates, and y coordinates as Eigen matrices
  std::tuple<MatXd, MatXd, MatXd> LoadTiff(std::filesystem::path path, std::array<double, 2> xlims,
                                           std::array<double, 2> ylims, double max_res) {
    GDALAllRegister();
    GDALDataset* dataset = (GDALDataset*)GDALOpen(path.string().c_str(), GA_ReadOnly);
    LUPNT_CHECK(dataset, "Failed to open TIFF", "LoadTiff");

    double gt[6];
    dataset->GetGeoTransform(gt);
    LUPNT_CHECK(gt, "Failed to get GeoTransform", "LoadTiff");

    int width_src = dataset->GetRasterXSize();
    int height_src = dataset->GetRasterYSize();
    double pixel_size_x = gt[1];
    double pixel_size_y = std::abs(gt[5]);

    double x_src_min = gt[0];
    double y_src_max = gt[3];
    double x_src_max = x_src_min + width_src * pixel_size_x;
    double y_src_min = y_src_max - height_src * pixel_size_y;

    double x_min = std::max(x_src_min, xlims[0]);
    double x_max = std::min(x_src_max, xlims[1]);
    double y_min = std::max(y_src_min, ylims[0]);
    double y_max = std::min(y_src_max, ylims[1]);

    int col_min = std::max(0, static_cast<int>(std::floor((x_min - gt[0]) / pixel_size_x)));
    int col_max = std::min(width_src, static_cast<int>(std::ceil((x_max - gt[0]) / pixel_size_x)));
    int row_min = std::max(0, static_cast<int>(std::floor((gt[3] - y_max) / pixel_size_y)));
    int row_max = std::min(height_src, static_cast<int>(std::ceil((gt[3] - y_min) / pixel_size_y)));

    int width = col_max - col_min;
    int height = row_max - row_min;

    int ds_x = std::max(1, static_cast<int>(std::ceil(max_res / pixel_size_x)));
    int ds_y = std::max(1, static_cast<int>(std::ceil(max_res / pixel_size_y)));

    int out_width = width / ds_x;
    int out_height = height / ds_y;

    std::vector<float> buffer(out_width * out_height);
    GDALRasterBand* band = dataset->GetRasterBand(1);

    GDALRasterIOExtraArg args;
    INIT_RASTERIO_EXTRA_ARG(args);
    args.eResampleAlg = GRIORA_NearestNeighbour;

    CPLErr err = band->RasterIO(GF_Read, col_min, row_min, width, height, buffer.data(), out_width,
                                out_height, GDT_Float32, 0, 0, &args);
    if (err != CE_None) LUPNT_CHECK(false, "Failed to read raster data", "LoadTiff");

    // Create Eigen matrices
    MatXd data(out_height, out_width);
    MatXd x_coords(out_height, out_width);
    MatXd y_coords(out_height, out_width);

    double origin_x = gt[0] + col_min * pixel_size_x;
    double origin_y = gt[3] - row_min * pixel_size_y;

    // Fill the matrices
    for (int r = 0; r < out_height; ++r) {
      for (int c = 0; c < out_width; ++c) {
        data(r, c) = static_cast<double>(buffer[r * out_width + c]);
        x_coords(r, c) = origin_x + (c + 0.5) * ds_x * pixel_size_x;
        y_coords(r, c) = origin_y - (r + 0.5) * ds_y * pixel_size_y;
      }
    }

    GDALClose(dataset);
    return {x_coords, y_coords, data};
  }
}  // namespace lupnt
