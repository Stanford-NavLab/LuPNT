#pragma once

#include <Eigen/Dense>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace lupnt {

  struct Crater {
    double lat;
    double lon;
    double diam_major;
    double diam_minor;
    double angle;
    std::string id;
  };

  class CraterDataLoader {
  public:
    static std::vector<Crater> LoadCraters(const std::string& filename
                                           = "lunar_crater_database_robbins_2018.csv",
                                           const std::pair<double, double>* latlims = nullptr,
                                           const std::pair<double, double>* longlims = nullptr,
                                           const std::pair<double, double>& diamlims = {0.0, 500.0},
                                           const double ellipse_limit = 1.5,
                                           const double arc_lims = 0.0);

    static void ExtractRobbinsDataset(const std::vector<Crater>& craters, Eigen::VectorXd& lat,
                                      Eigen::VectorXd& lon, Eigen::VectorXd& major,
                                      Eigen::VectorXd& minor, Eigen::VectorXd& psi,
                                      std::vector<std::string>& crater_id, bool radians = true);
  };

}  // namespace lupnt
