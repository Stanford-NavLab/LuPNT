#pragma once

#include <yaml-cpp/yaml.h>

#include "lupnt/core/definitions.h"

namespace YAML {
  /// @brief `yaml-cpp` `convert<>` specialization enabling `node.as<lupnt::Real>()` /
  /// `node = some_real` for LuPNT's autodiff scalar type.
  ///
  /// Lets simulation YAML config files (loaded via `lupnt::GetYamlConfig`,
  /// `lupnt/core/file.h`) be read directly into `Real` parameters -- e.g.
  /// dynamics/device config values such as `config["dt"].as<Real>()`,
  /// `config["mass"].as<Real>()` (numerical_orbit_dynamics.cc),
  /// `cfg["lat_ref"].as<Real>()` (rover.cc) -- without callers needing to
  /// round-trip through `double` themselves.
  template <> struct convert<lupnt::Real> {
    /// @brief Encode a `Real` as a YAML scalar node (its underlying `double` value).
    static Node encode(const lupnt::Real& rhs) {
      Node node;
      node = static_cast<double>(rhs);
      return node;
    }

    /// @brief Decode a YAML scalar node into a `Real`. Returns false (decode
    /// failure) if `node` is not a scalar.
    static bool decode(const Node& node, lupnt::Real& rhs) {
      if (!node.IsScalar()) return false;
      rhs = node.as<double>();
      return true;
    }
  };

  /// @brief `yaml-cpp` `convert<>` specialization enabling
  /// `node.as<Eigen::Matrix<Scalar, Rows, Cols>>()` for fixed- or
  /// dynamic-size Eigen vectors/matrices (e.g. `Vec3`, `Vec6`, `Mat3`,
  /// `VecXd`).
  ///
  /// Lets simulation YAML configs specify vector/matrix-valued parameters
  /// (e.g. initial position/velocity, attitude matrices, covariance/noise
  /// matrices) as nested sequences -- `[x, y, z]` for a vector or
  /// `[[..],[..],[..]]` for a matrix -- which are parsed directly into the
  /// corresponding Eigen type via `node.as<Vec3>()` etc. Row/column vectors
  /// (`Rows==1` or `Cols==1`) are read from a flat sequence; otherwise `node`
  /// must be a sequence of equal-length row sequences.
  template <typename Scalar, int Rows, int Cols> struct convert<Eigen::Matrix<Scalar, Rows, Cols>> {
    /// @brief Encode an Eigen matrix as a YAML sequence of row sequences.
    static Node encode(const Eigen::Matrix<Scalar, Rows, Cols>& mat) {
      Node node;
      for (int i = 0; i < mat.rows(); ++i) {
        Node row;
        for (int j = 0; j < mat.cols(); ++j) row.push_back(mat(i, j));
        node.push_back(row);
      }
      return node;
    }

    /// @brief Decode a YAML sequence (flat for vectors, nested for matrices)
    /// into `mat`, resizing dynamic dimensions as needed. Returns false on
    /// shape mismatch or if `node` is not a sequence.
    static bool decode(const Node& node, Eigen::Matrix<Scalar, Rows, Cols>& mat) {
      if (!node.IsSequence()) return false;

      // Handle vector, row or column, as a flat sequence.
      constexpr bool is_row_vector = (Rows == 1 && Cols != 1);
      constexpr bool is_col_vector = (Cols == 1 && Rows != 1);
      constexpr bool is_vector = (is_row_vector || is_col_vector);

      if (is_vector) {
        int vec_size = static_cast<int>(node.size());

        if (Rows == 1 && Cols == Eigen::Dynamic) {
          mat.derived().resize(1, vec_size);
        } else if (Cols == 1 && Rows == Eigen::Dynamic) {
          mat.derived().resize(vec_size, 1);
        } else if (Rows == 1 && Cols != Eigen::Dynamic && Cols != vec_size) {
          return false;
        } else if (Cols == 1 && Rows != Eigen::Dynamic && Rows != vec_size) {
          return false;
        }

        for (int i = 0; i < vec_size; ++i) {
          if (is_row_vector) {
            mat(0, i) = node[i].as<Scalar>();
          } else {
            mat(i, 0) = node[i].as<Scalar>();
          }
        }
        return true;
      }

      // Matrix case.
      if (Rows != Eigen::Dynamic && node.size() != static_cast<size_t>(Rows)) {
        return false;
      }

      int n_rows = static_cast<int>(node.size());
      int n_cols = 0;

      if (n_rows > 0) {
        if (!node[0].IsSequence()) return false;
        n_cols = static_cast<int>(node[0].size());
        if (Cols != Eigen::Dynamic && n_cols != Cols) return false;
      }

      if (Rows == Eigen::Dynamic && Cols == Eigen::Dynamic) {
        mat.derived().resize(n_rows, n_cols);
      } else if (Rows == Eigen::Dynamic) {
        mat.derived().resize(n_rows, mat.cols());
      } else if (Cols == Eigen::Dynamic) {
        mat.derived().resize(mat.rows(), n_cols);
      }

      for (int i = 0; i < n_rows; ++i) {
        const Node& row = node[i];
        if (!row.IsSequence()) return false;
        if (Cols != Eigen::Dynamic && static_cast<int>(row.size()) != Cols) return false;

        for (int j = 0; j < static_cast<int>(row.size()); ++j) {
          mat(i, j) = row[j].as<Scalar>();
        }
      }

      return true;
    }
  };
}  // namespace YAML
