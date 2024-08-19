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
  public:
    File(const std::string& filepath);
    ~File();
    template <typename T> File& operator<<(const T& data);

  private:
    std::ofstream file;
  };

  template <typename T> class Timestamped {
  public:
    Timestamped(double timestamp, const T& data);
    double GetTimestamp() const;
    const T& GetData() const;

  private:
    double timestamp;
    T data;
  };

  class DataHistory {
  public:
    template <typename VecType>
    void AddData(const std::string& key, double timestamp, const VecType& data);
    const std::vector<Timestamped<VecXd>>& GetData(const std::string& key) const;
    void AddHeader(const std::string& key, const std::string& header);
    const std::map<std::string, std::vector<Timestamped<VecXd>>>& GetData() const;
    const std::map<std::string, std::vector<std::string>>& GetHeaders() const;

  private:
    std::map<std::string, std::vector<Timestamped<VecXd>>> historyData;
    std::map<std::string, std::vector<std::string>> headers;
  };

  class FileWriter {
  public:
    FileWriter(const std::filesystem::path& basePath, const bool make_dirs = false);
    void WriteData(const DataHistory& dataHistory);

  private:
    std::filesystem::path basePath;
    Eigen::IOFormat fmt;
  };

  size_t CountLines(const std::filesystem::path& filepath);

}  // namespace lupnt
