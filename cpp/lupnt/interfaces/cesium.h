#pragma once

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
#include <filesystem>
#include <mutex>
#include <nlohmann/json.hpp>
#include <string>
#include <thread>
#include <vector>

#include "lupnt/core/constants.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  struct Entity {
    std::string id;
    std::string name;
    VecX times;
    MatX3 positions;
    std::vector<int> color;
    std::string description;
    BodyId body_id;
    std::string initial_time_utc;
    int size;
  };

  class CesiumViewer {
  public:
    /// @brief Construct a Cesium-based 3D visualization server.
    ///
    /// Created by `Simulation::Init` (see `lupnt/simulations/simulation.h`,
    /// `GetCesiumViewer`) when the simulation config enables Cesium output;
    /// agents (e.g. `Satellite`) then call `AddEntity` each step to stream
    /// their trajectory to the served CZML viewer. Requires the
    /// `CESIUM_TOKEN` environment variable to be set; copies and patches the
    /// `CesiumViewer/index.html` template (with the token substituted) into
    /// the local cache directory.
    ///
    /// @param port Local TCP port to serve the viewer on
    /// @param host Host/IP address to bind the server to
    CesiumViewer(int port = 8080, const std::string& host = "127.0.0.1");

    /// @brief Stops the server (if running) and releases resources.
    ~CesiumViewer();

    /// @brief Start the Cesium web server on a background thread.
    ///
    /// Called once after construction (typically by `Simulation::Init`) so
    /// the viewer is reachable at `http://<host>:<port>` for the remainder of
    /// the simulation run.
    void Run();

    /// @brief Register (or replace) a trajectory entity to be displayed in
    /// the Cesium viewer.
    ///
    /// Called once per agent (and again whenever its propagated trajectory is
    /// updated) to push the agent's time history of positions for rendering
    /// as a CZML path; served on request via `/entities/<id>.czml`.
    ///
    /// @param times       Sample times [s, elapsed since the simulation's
    ///                    initial UTC epoch]
    /// @param positions   Body-fixed (`body_id`'s "FIXED" frame) positions
    ///                    [m], one row per sample (x, y, z columns)
    /// @param entity_id   Unique identifier for the entity (used in CZML and
    ///                    the `/entities/<id>.czml` route)
    /// @param name        Display name shown in the Cesium viewer
    /// @param color       RGB color `{r, g, b}` (0-255) used for the point
    ///                     and path
    /// @param description Text description shown in the entity's info box
    /// @param body_id     Central body the positions are expressed relative
    ///                    to (sets the CZML reference body)
    /// @param size        Point marker size [pixels]
    void AddEntity(const VecX& times, const MatX3& positions,
                   const std::string& entity_id = "Entity", const std::string& name = "Entity",
                   const std::vector<int>& color = {255, 0, 255},
                   const std::string& description = "Entity orbit", BodyId body_id = BodyId::MOON,
                   int size = 5);

    /// @brief Remove all registered entities from the viewer.
    void ClearEntities();

    /// @brief Stop the background server thread started by `Run`.
    void Stop();

  private:
    crow::SimpleApp app_;
    std::filesystem::path static_path_;
    std::string host_;
    int port_;
    std::vector<Entity> entities_;
    std::mutex mutex_;
    std::thread server_thread_;
    bool running_ = false;

    /// @brief Register the HTTP routes served by `app_` (index page, entity
    /// list, time range, and per-entity CZML).
    void SetupRoutes();

    /// @brief Build the CZML JSON document (as a string) describing
    /// `entity`'s point, label, and path, for the `/entities/<id>.czml` route.
    std::string CreateCzmlDataForEntity(const Entity& entity);

    /// @brief Parse an ISO 8601 UTC timestamp (`YYYY-MM-DDTHH:MM:SSZ`) into a
    /// `std::chrono::system_clock::time_point`.
    static std::chrono::system_clock::time_point ParseIso8601(const std::string& s);
  };

}  // namespace lupnt
