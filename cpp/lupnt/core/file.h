/**
 * @file file.h
 * @author Stanford NAV LAB
 * @brief File access utils
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <highfive/H5Easy.hpp>
#include <iostream>
#include <map>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

#include "lupnt/core/definitions.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {

  std::filesystem::path GetBaseDir();
  std::filesystem::path GetDataPath();
  std::filesystem::path GetInputPath(const std::string& input_file);
  YAML::Node GetYamlConfig(const std::filesystem::path& filepath);

  std::filesystem::path GetOutputDir(const std::string& dirname);
  std::filesystem::path SetOutputDir(const std::string& dirname);
  H5Easy::File GetH5File(std::filesystem::path cache_path, bool recompute);

  template <typename Derived> void Dump(H5Easy::File& file, const std::string& name,
                                        const Eigen::DenseBase<Derived>& matrix,
                                        H5Easy::DumpMode mode = H5Easy::DumpMode::Overwrite) {
    using Scalar = typename Derived::Scalar;
    if constexpr (std::is_same_v<Scalar, bool>) {
      std::vector<bool> bool_data(matrix.size());
      for (int i = 0; i < matrix.size(); ++i) bool_data[i] = matrix(i);
      H5Easy::dump(file, name, bool_data, mode);
    } else {
      MatXd tmp = matrix.template cast<double>();
      H5Easy::dump(file, name, tmp, mode);
    }
  }

  void Dump(H5Easy::File& file, const std::string& name, Real value,
            H5Easy::DumpMode mode = H5Easy::DumpMode::Overwrite);

  template <int N, int M, typename T, typename Func, typename... Args>
  Matrix<T, N, M> LoadOrRecompute(const std::string& name, H5Easy::File& cache_file, bool recompute,
                                  Func func, Args&&... args) {
    Matrix<T, N, M> result;
    bool loaded = false;
#pragma omp critical
    {
      if (!recompute && cache_file.exist(name)) {
        if constexpr (std::is_same_v<T, bool>) {
          std::vector<uint8_t> bool_data = H5Easy::load<std::vector<uint8_t>>(cache_file, name);
          result = Eigen::Map<const Matrix<bool, N, M>>(
              reinterpret_cast<const bool*>(bool_data.data()), N, M);
        } else {
          result = H5Easy::load<Matrix<double, N, M>>(cache_file, name);
        }
        loaded = true;
      }
    }
    if (loaded) return result.template cast<T>();

    result = func(std::forward<Args>(args)...).template cast<T>();
#pragma omp critical
    {
      if constexpr (std::is_same_v<T, bool>) {
        std::vector<uint8_t> bool_data(result.size());
        for (int i = 0; i < result.size(); ++i) {
          bool_data[i] = static_cast<uint8_t>(result.data()[i]);
        }
        H5Easy::dump(cache_file, name, bool_data, H5Easy::DumpMode::Overwrite);
      } else {
        Dump(cache_file, name, result, H5Easy::DumpMode::Overwrite);
      }
    }
    return result;
  }

  std::optional<std::filesystem::path> FindFileInDir(const std::filesystem::path& base_path,
                                                     std::string_view filename);
  std::filesystem::path GetFilePath(std::string_view filename);
  std::filesystem::path GetCspiceKernelDir();
  std::filesystem::path GetAsciiKernelDir();
  std::filesystem::path GetConfigFileDir();

  std::chrono::time_point<std::chrono::high_resolution_clock> GetSystemTime();
  std::string PrintDuration(const std::chrono::duration<double>& duration, int precision = 3);

  template <typename T> T OpenFile(const std::filesystem::path& filepath);
  size_t CountLines(const std::filesystem::path& filepath);

}  // namespace lupnt
