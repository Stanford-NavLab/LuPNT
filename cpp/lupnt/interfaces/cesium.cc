#include "lupnt/interfaces/cesium.h"

// Suppress warnings from crow third-party library
#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wuninitialized"
#  pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wuninitialized"
#endif

#include <crow.h>

#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif
#ifdef __clang__
#  pragma clang diagnostic pop
#endif

#include <Eigen/Dense>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <nlohmann/json.hpp>
#include <optional>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "lupnt/conversions/time_conversions.h"
#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"
#include "lupnt/core/error.h"
#include "lupnt/core/file.h"
#include "lupnt/core/logger.h"

using json = nlohmann::json;

namespace lupnt {
  CesiumViewer::CesiumViewer(int port, const std::string& host) : host_(host), port_(port) {
    static_path_ = GetBaseDir() / ".cache" / "cesium";
    std::filesystem::create_directories(static_path_);

    const char* cesium_token = std::getenv("CESIUM_TOKEN");
    if (!cesium_token)
      LUPNT_CHECK(false, "CESIUM_TOKEN environment variable is not set", "CesiumViewer");

    // Copy index.html template and replace token
    std::filesystem::path index_template = GetDataPath() / "CesiumViewer" / "index.html";
    std::ifstream in(index_template);
    std::string index_content((std::istreambuf_iterator<char>(in)),
                              std::istreambuf_iterator<char>());
    size_t pos = index_content.find("{{ CESIUM_TOKEN }}");
    if (pos != std::string::npos) {
      index_content.replace(pos, std::string("{{ CESIUM_TOKEN }}").length(), cesium_token);
    } else {
      LUPNT_CHECK(false, "Failed to replace CESIUM_TOKEN in index.html", "CesiumViewer");
    }
    std::ofstream out(static_path_ / "index.html");
    out << index_content;
    out.close();
    Logger::Debug("Copied index.html to " + static_path_.string(), "CesiumViewer");

    app_.loglevel(crow::LogLevel::Warning);
    SetupRoutes();
  }

  CesiumViewer::~CesiumViewer() { Stop(); }

  void CesiumViewer::Stop() {
    if (!running_) return;
    app_.stop();
    Logger::Info("Stopping CesiumViewer", "CesiumViewer");
    if (server_thread_.joinable()) server_thread_.join();
    running_ = false;
  }

  void CesiumViewer::Run() {
    Logger::Info("Running on http://" + host_ + ":" + std::to_string(port_), "CesiumViewer");
    running_ = true;
    server_thread_ = std::thread([this]() { app_.port(port_).multithreaded().run(); });
  }

  void CesiumViewer::AddEntity(const VecX& times, const MatX3& positions,
                               const std::string& entity_id, const std::string& name,
                               const std::vector<int>& color, const std::string& description,
                               BodyId body_id, int size) {
    std::string frame = "FIXED";
    LUPNT_CHECK(times.size() == positions.rows(),
                "Times and positions arrays must have the same length", "CesiumViewer");
    std::string initial_time_utc = TimeToGregorianString(UtcToTai(GetLupntEpoch()));

    Entity entity{entity_id,        name, times, positions, color, description, body_id,
                  initial_time_utc, size};

    std::lock_guard<std::mutex> lock(mutex_);
    entities_.push_back(entity);
    Logger::Info("Added entity " + entity_id, "CesiumViewer");
  }

  void CesiumViewer::ClearEntities() {
    std::lock_guard<std::mutex> lock(mutex_);
    entities_.clear();
    Logger::Info("Cleared entities", "CesiumViewer");
  }

  void CesiumViewer::SetupRoutes() {
    Logger::Debug("Setting up routes", "CesiumViewer");

    CROW_ROUTE(app_, "/")([this]() {
      std::ifstream file(static_path_ / "index.html");
      if (!file) return crow::response(404);
      std::ostringstream contents;
      contents << file.rdbuf();
      Logger::Debug("Serving index.html", "CesiumViewer");
      return crow::response(contents.str());
    });

    CROW_ROUTE(app_, "/entities/all")([this]() {
      json entities_info = json::array();
      std::lock_guard<std::mutex> lock(mutex_);
      for (const auto& entity : entities_) {
        entities_info.push_back({{"entity_id", entity.id}, {"body_id", entity.body_id}});
      }
      json result;
      result["entities"] = entities_info;
      Logger::Debug("Serving entities", "CesiumViewer");
      return crow::response(result.dump());
    });

    CROW_ROUTE(app_, "/entities/timerange")([this]() {
      std::lock_guard<std::mutex> lock(mutex_);
      if (entities_.empty()) {
        return crow::response(R"({"start_time": null, "end_time": null})");
      }
      std::optional<std::chrono::system_clock::time_point> global_start, global_end;
      for (const auto& entity : entities_) {
        auto start = ParseIso8601(entity.initial_time_utc);
        double last_time = (entity.times.size() == 0)
                               ? 0.0
                               : static_cast<double>(entity.times(entity.times.size() - 1));
        auto stop = start
                    + std::chrono::duration_cast<std::chrono::system_clock::duration>(
                        std::chrono::duration<double>(last_time));
        if (!global_start || start < *global_start) global_start = start;
        if (!global_end || stop > *global_end) global_end = stop;
      }
      auto to_iso = [](const std::chrono::system_clock::time_point& tp) {
        std::time_t t = std::chrono::system_clock::to_time_t(tp);
        std::tm tm = *std::gmtime(&t);
        char buf[32];
        std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
        return std::string(buf) + "Z";
      };
      json result;
      result["start_time"] = to_iso(*global_start);
      result["end_time"] = to_iso(*global_end);
      Logger::Debug("Serving timerange", "CesiumViewer");
      return crow::response(result.dump());
    });

    CROW_ROUTE(app_, "/entities/<path>")
    ([this](const crow::request&, crow::response& res, std::string entity_path) {
      Logger::Debug("Serving entity " + entity_path, "CesiumViewer");
      if (entity_path.size() > 5 && entity_path.substr(entity_path.size() - 5) == ".czml") {
        std::string entity_id = entity_path.substr(0, entity_path.size() - 5);
        std::optional<Entity> entity_opt;
        {
          std::lock_guard<std::mutex> lock(mutex_);
          for (const auto& e : entities_) {
            if (e.id == entity_id) {
              entity_opt = e;
              break;
            }
          }
        }
        if (!entity_opt) {
          res.code = 404;
          res.write("Entity not found");
          res.end();
          return;
        }
        std::string czml = CreateCzmlDataForEntity(*entity_opt);
        res.set_header("Content-Type", "application/json");
        res.write(czml);
        res.end();
      } else {
        res.code = 404;
        res.end();
      }
    });

    CROW_ROUTE(app_, "/<path>")
    ([](const crow::request&, crow::response& res, std::string path) {
      Logger::Debug("Serving path " + path, "CesiumViewer");
      res.code = 404;
      res.end();
    });
  }

  std::string CesiumViewer::CreateCzmlDataForEntity(const Entity& entity) {
    auto start = ParseIso8601(entity.initial_time_utc);
    double last_time = (entity.times.size() == 0)
                           ? 0.0
                           : static_cast<double>(entity.times(entity.times.size() - 1));
    auto stop = start
                + std::chrono::duration_cast<std::chrono::system_clock::duration>(
                    std::chrono::duration<double>(last_time));
    auto to_iso = [](const std::chrono::system_clock::time_point& tp) {
      std::time_t t = std::chrono::system_clock::to_time_t(tp);
      std::tm tm = *std::gmtime(&t);
      char buf[32];
      std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
      return std::string(buf) + "Z";
    };

    std::vector<double> flat_positions;
    std::vector<std::string> time_strings;
    for (Eigen::Index i = 0; i < entity.times.size(); ++i) {
      auto t = start
               + std::chrono::duration_cast<std::chrono::system_clock::duration>(
                   std::chrono::duration<double>(entity.times(i)));
      time_strings.push_back(to_iso(t));
      flat_positions.push_back(entity.positions(i, 0) * M_KM);
      flat_positions.push_back(entity.positions(i, 1) * M_KM);
      flat_positions.push_back(entity.positions(i, 2) * M_KM);
    }

    // Compose position data as [time, x, y, z, ...]
    std::vector<json> position_data;
    for (size_t i = 0; i < time_strings.size(); ++i) {
      position_data.push_back(time_strings[i]);
      position_data.push_back(flat_positions[i * 3 + 0]);
      position_data.push_back(flat_positions[i * 3 + 1]);
      position_data.push_back(flat_positions[i * 3 + 2]);
    }

    // Use color with fallback to magenta if not provided
    std::vector<int> color = entity.color;
    if (color.size() != 3) color = {255, 0, 255};

    json czml_entity = {
        {"id", entity.id},
        {"name", entity.name},
        {"description", entity.description},
        {"availability", to_iso(start) + "/" + to_iso(stop)},
        {"label",
         {
             {"fillColor", {{"rgba", {255, 255, 255, 255}}}},
             {"font", "13pt Lucida Console"},
             {"horizontalOrigin", "LEFT"},
             {"outlineColor", {{"rgba", {0, 0, 0, 255}}}},
             {"outlineWidth", 3},
             {"pixelOffset", {{"cartesian2", {20, 0}}}},
             {"style", "FILL_AND_OUTLINE"},
             {"text", entity.name},
         }},
        {"position",
         {
             {"cartesian", position_data},
             {"interpolationAlgorithm", "LAGRANGE"},
             {"interpolationDegree", 1},
             {"referenceFrame", "FIXED"},
         }},
        {"point",
         {
             {"pixelSize", entity.size},
             {"color", {{"rgba", {color[0], color[1], color[2], 255}}}},
             {"outlineColor", {{"rgba", {color[0], color[1], color[2], 255}}}},
             {"outlineWidth", 2},
             {"show", true},
         }},
    };

    if (entity.times.size() > 1) {
      czml_entity["path"] = {
          {"show", true},
          {"width", 2},
          {"material",
           {{"solidColor", {{"color", {{"rgba", {color[0], color[1], color[2], 128}}}}}}}},
          {"resolution", 120},
      };
    }

    json czml_doc = {
        {{"id", "document"}, {"name", "Satellite Visualization"}, {"version", "1.0"}},
        czml_entity,
    };
    return czml_doc.dump(2);
  }

  // Parse ISO8601 string (YYYY-MM-DDTHH:MM:SSZ) to time_point
  std::chrono::system_clock::time_point CesiumViewer::ParseIso8601(const std::string& s) {
    std::tm tm = {};
    std::istringstream ss(s);
    ss >> std::get_time(&tm, "%Y-%m-%dT%H:%M:%S");
    if (s.back() == 'Z') {
      // UTC
    }
    return std::chrono::system_clock::from_time_t(std::mktime(&tm));
  }

}  // namespace lupnt
