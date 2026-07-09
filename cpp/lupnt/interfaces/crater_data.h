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
    /// @brief Load and filter lunar crater records from a Robbins-format crater database CSV.
    ///
    /// Used by surface/rover environment models to build a catalog of lunar craters as
    /// landmark/feature candidates (e.g. for terrain-relative navigation or crater-detection
    /// measurement models). Reads `filename` from the bundled surface data directory and keeps
    /// only craters that pass all of the supplied region/size/shape filters.
    ///
    /// @param filename     Name of the crater database CSV file (Robbins 2018 format), resolved
    ///                      relative to the bundled lunar surface data directory.
    /// @param latlims      Optional [min, max] latitude bounds [deg]; craters outside this range
    ///                      are excluded. If nullptr, no latitude filtering is applied.
    /// @param longlims     Optional [min, max] longitude bounds [deg]; craters outside this range
    ///                      are excluded. If nullptr, no longitude filtering is applied.
    /// @param diamlims     [min, max] crater diameter bounds [km]; craters with major-axis
    ///                      diameter outside this range are excluded (default: [0, 500] km).
    /// @param ellipse_limit Maximum allowed major/minor diameter ratio; craters more elongated
    ///                      than this (i.e. less circular) are excluded (default: 1.5).
    /// @param arc_lims     Minimum required fraction of the crater rim arc that was imaged/
    ///                      digitized in the source database; craters below this threshold are
    ///                      excluded (default: 0.0, i.e. no filtering).
    /// @return             Filtered list of Crater records (lat/lon [deg], major/minor diameters
    ///                      [km], ellipse angle, and database ID).
    static std::vector<Crater> LoadCraters(const std::string& filename
                                           = "lunar_crater_database_robbins_2018.csv",
                                           const std::pair<double, double>* latlims = nullptr,
                                           const std::pair<double, double>* longlims = nullptr,
                                           const std::pair<double, double>& diamlims = {0.0, 500.0},
                                           const double ellipse_limit = 1.5,
                                           const double arc_lims = 0.0);

    /// @brief Split a list of Crater records into parallel column vectors for vectorized use.
    ///
    /// Used to convert the per-crater records returned by LoadCraters into the flat Eigen vector
    /// form expected by crater-matching/landmark-based measurement and navigation routines.
    ///
    /// @param craters   Crater records to extract, as returned by LoadCraters.
    /// @param lat       Output latitude of each crater [deg, or rad if `radians` is true].
    /// @param lon       Output longitude of each crater [deg, or rad if `radians` is true].
    /// @param major     Output major-axis diameter of each crater [km].
    /// @param minor     Output minor-axis diameter of each crater [km].
    /// @param psi       Output ellipse orientation angle of each crater [deg, or rad if
    ///                  `radians` is true].
    /// @param crater_id Output database ID string of each crater.
    /// @param radians   If true, `lat`/`lon`/`psi` are converted to radians; otherwise left in
    ///                  degrees (default: true).
    static void ExtractRobbinsDataset(const std::vector<Crater>& craters, Eigen::VectorXd& lat,
                                      Eigen::VectorXd& lon, Eigen::VectorXd& major,
                                      Eigen::VectorXd& minor, Eigen::VectorXd& psi,
                                      std::vector<std::string>& crater_id, bool radians = true);
  };

}  // namespace lupnt
