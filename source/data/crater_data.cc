#include "lupnt/data/crater_data.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "lupnt/core/user_file_path.h"
#include "lupnt/numerics/string_utils.h"

namespace lupnt {

  // Function to load craters based on various limits
  std::vector<Crater> CraterDataLoader::LoadCraters(const std::string& filename,
                                                    const std::pair<double, double>* latlims,
                                                    const std::pair<double, double>* longlims,
                                                    const std::pair<double, double>& diamlims,
                                                    const double ellipse_limit,
                                                    const double arc_lims) {
    auto path = GetFilePath(filename);
    auto data = ReadCSV(path);
    std::vector<Crater> craters;

    for (const auto& row : data) {
      if (row.empty()) continue;

      try {
        double arc_img = std::stod(row.at(10));
        if (arc_img <= arc_lims) continue;

        double lat = std::stod(row.at(1));
        double lon = std::stod(row.at(2));
        double diam_major = std::stod(row.at(4));
        double diam_minor = std::stod(row.at(5));
        double diam_circ = std::stod(row.at(3));

        if (latlims && (lat < latlims->first || lat > latlims->second)) continue;
        if (longlims && (lon < longlims->first || lon > longlims->second)) continue;
        if (diam_circ < diamlims.first || diam_circ > diamlims.second) continue;

        double ellipticity = diam_major / diam_minor;
        if (ellipticity > ellipse_limit) continue;

        craters.push_back({lat, lon, diam_major, diam_minor, std::stod(row.at(6)), row.at(0)});
      } catch (const std::out_of_range& e) {
        // Handle case where a field is missing
        continue;
      } catch (const std::invalid_argument& e) {
        // Handle case where conversion failed (empty or non-numeric field)
        continue;
      }
    }

    return craters;
  }

  // Function to extract dataset
  void CraterDataLoader::ExtractRobbinsDataset(const std::vector<Crater>& craters,
                                               Eigen::VectorXd& lat, Eigen::VectorXd& lon,
                                               Eigen::VectorXd& major, Eigen::VectorXd& minor,
                                               Eigen::VectorXd& psi,
                                               std::vector<std::string>& crater_id, bool radians) {
    size_t size = craters.size();
    lat.resize(size);
    lon.resize(size);
    major.resize(size);
    minor.resize(size);
    psi.resize(size);
    crater_id.resize(size);

    for (size_t i = 0; i < size; ++i) {
      lat[i] = radians ? craters[i].lat * M_PI / 180.0 : craters[i].lat;
      lon[i] = radians ? craters[i].lon * M_PI / 180.0 : craters[i].lon;
      major[i] = craters[i].diam_major;
      minor[i] = craters[i].diam_minor;
      psi[i] = radians ? craters[i].angle * M_PI / 180.0 : craters[i].angle;
      crater_id[i] = craters[i].id;
    }
  }

}  // namespace lupnt