#include "lupnt/core/file.h"

namespace lupnt {

  File::File(const std::string& filepath) {
    file.open(filepath);
    if (!file.is_open()) {
      throw std::runtime_error("Could not open file: " + filepath);
    }
  }

  File::~File() {
    if (file.is_open()) {
      file.close();
    }
  }

  template <typename T> File& File::operator<<(const T& data) {
    file << data;
    return *this;
  }

  template <typename T> Timestamped<T>::Timestamped(double timestamp, const T& data)
      : timestamp(timestamp), data(data) {}

  template <typename T> double Timestamped<T>::GetTimestamp() const { return timestamp; }
  template <typename T> const T& Timestamped<T>::GetData() const { return data; }

  template <typename VecType>
  void DataHistory::AddData(const std::string& key, double timestamp, const VecType& data) {
    if (historyData.find(key) == historyData.end()) {
      historyData[key] = {};
    }
    historyData[key].push_back(Timestamped<VecXd>(timestamp, data.template cast<double>()));
  }

  const std::vector<Timestamped<VecXd>>& DataHistory::GetData(const std::string& key) const {
    auto it = historyData.find(key);
    if (it != historyData.end()) {
      return it->second;
    }
    throw std::runtime_error("No data found for the given key.");
  }

  void DataHistory::AddHeader(const std::string& key, const std::string& header) {
    headers[key].push_back(header);
  }

  const std::map<std::string, std::vector<Timestamped<VecXd>>>& DataHistory::GetData() const {
    return historyData;
  }

  const std::map<std::string, std::vector<std::string>>& DataHistory::GetHeaders() const {
    return headers;
  }

  FileWriter::FileWriter(const std::filesystem::path& basePath, const bool make_dirs)
      : basePath(basePath) {
    if (make_dirs) {
      std::filesystem::create_directories(basePath);
    }
    if (!std::filesystem::exists(basePath)) {
      throw std::runtime_error("Base path does not exist: " + basePath.string());
    }
    for (const auto& entry : std::filesystem::directory_iterator(basePath)) {
      if (entry.path().extension() == ".csv") {
        std::filesystem::remove(entry.path());
      }
    }
  }

  void FileWriter::WriteData(const DataHistory& dataHistory) {
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

  size_t CountLines(const std::filesystem::path& filepath) {
    std::ifstream file(filepath);
    assert(file.is_open() && "Unable to open file");

    size_t lineCount = 0;
    std::string line;
    while (std::getline(file, line)) {
      ++lineCount;
    }
    file.close();
    return lineCount;
  }

}  // namespace lupnt
