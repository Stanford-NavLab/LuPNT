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

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "lupnt/numerics/math_utils.h"

namespace lupnt {

class File {
 private:
  std::ofstream file;

 public:
  File(const std::string& filepath) {
    file.open(filepath);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open file: " + filepath);
    }
  }

  ~File() {
    if (file.is_open()) {
      file.close();
    }
  }

  // Operator <<
  template <typename T>
  File& operator<<(const T& data) {
    file << data;
    return *this;
  }
};

template <typename T>
class Timestamped {
 public:
  Timestamped(double timestamp, const T& data)
      : timestamp(timestamp), data(data) {}

  double GetTimestamp() const { return timestamp; }
  const T& GetData() const { return data; }

 private:
  double timestamp;
  T data;
};

class DataHistory {
 public:
  template <typename VectorType>
  void AddData(const std::string& key, double timestamp,
               const VectorType& data) {
    if (historyData.find(key) == historyData.end()) {
      historyData[key] = {};
    }
    historyData[key].push_back(
        Timestamped<VectorXd>(timestamp, data.template cast<double>()));
  }

  const std::vector<Timestamped<VectorXd>>& GetData(
      const std::string& key) const {
    auto it = historyData.find(key);
    if (it != historyData.end()) {
      return it->second;
    }
    throw std::runtime_error("No data found for the given key.");
  }

  void AddHeader(const std::string& key, const std::string& header) {
    headers[key].push_back(header);
  }

  const std::map<std::string, std::vector<Timestamped<VectorXd>>>& GetData()
      const {
    return historyData;
  }

  const std::map<std::string, std::vector<std::string>>& GetHeaders() const {
    return headers;
  }

 private:
  std::map<std::string, std::vector<Timestamped<VectorXd>>> historyData;
  std::map<std::string, std::vector<std::string>> headers;
};

class FileWriter {
 public:
  FileWriter(const std::filesystem::path& basePath,
             const bool make_dirs = false)
      : basePath(basePath) {
    // Create the base path if it does not exist
    if (make_dirs) {
      std::filesystem::create_directories(basePath);
    }
    // Check if the base path exists
    if (!std::filesystem::exists(basePath)) {
      throw std::runtime_error("Base path does not exist: " +
                               basePath.string());
    }
    // Remove all csv files in the directory
    for (const auto& entry : std::filesystem::directory_iterator(basePath)) {
      if (entry.path().extension() == ".csv") {
        std::filesystem::remove(entry.path());
      }
    }
  }

  void WriteData(const DataHistory& dataHistory) {
    for (const auto& [key, data] : dataHistory.GetData()) {
      std::filesystem::path filepath = basePath / (key + ".csv");

      File file(filepath.string());
      // Check headers
      auto it = dataHistory.GetHeaders().find(key);
      if (it != dataHistory.GetHeaders().end()) {
        file << "t";
        for (const auto& header : it->second) {
          file << "," << header;
        }
        file << "\n";
      }
      for (const auto& timestampedData : data) {
        file << timestampedData.GetTimestamp();
        if (timestampedData.GetData().size() > 0) {
          file << "," << timestampedData.GetData().transpose().format(fmt);
        }
        file << "\n";
      }
    }
  }

 private:
  std::filesystem::path basePath;
  Eigen::IOFormat fmt{Eigen::FullPrecision, Eigen::DontAlignCols, ",", "\n"};
};

}  // namespace lupnt