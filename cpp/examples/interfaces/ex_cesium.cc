// Builds a minimal CZML packet for orbit visualization in Cesium. Most of the
// server code is left commented as a template for local web visualization.
#include <lupnt/lupnt.h>

using namespace lupnt;
using namespace matplot;

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <crow.h>

json CreateCzml(const MatXd &rvs, const VecXd &tspan) {
  // CZML uses a document packet followed by one or more entity packets. The
  // position array stores [seconds since epoch, x, y, z] tuples.
  json czml = json::array();

  json document;
  document["id"] = "document";
  document["version"] = "1.0";
  czml.push_back(document);

  json entity;
  entity["id"] = "Satellite/Sat";
  entity["name"] = "Sat";
  entity["availability"] = "2012-03-15T10:00:00Z/2012-03-16T10:00:00Z";
  entity["description"] = "Sat orbit";

  json position;
  position["interpolationAlgorithm"] = "LAGRANGE";
  position["interpolationDegree"] = 5;
  position["referenceFrame"] = "INERTIAL";
  position["epoch"] = "2012-03-15T10:00:00Z";

  json positions = json::array();
  for (size_t i = 0; i < tspan.size(); ++i) {
    positions.push_back(tspan[i] - tspan[0]);
    positions.push_back(rvs(i, 0) * 1000);
    positions.push_back(rvs(i, 1) * 1000);
    positions.push_back(rvs(i, 2) * 1000);
  }

  position["cartesian"] = positions;
  entity["position"] = position;

  entity["label"] = {{"fillColor", {{"rgba", {255, 0, 255, 255}}}},
                     {"font", "11pt Lucida Console"},
                     {"horizontalOrigin", "LEFT"},
                     {"verticalOrigin", "CENTER"},
                     {"outlineColor", {{"rgba", {0, 0, 0, 255}}}},
                     {"outlineWidth", 2},
                     {"pixelOffset", {{"cartesian2", {12, 0}}}},
                     {"style", "FILL_AND_OUTLINE"},
                     {"text", "ISS"},
                     {"show", true}};

  entity["path"] = {{"show", true},
                    {"width", 1},
                    {"material", {"solidColor", {"color", {"rgba", {255, 0, 255, 255}}}}},
                    {"resolution", 120}};

  czml.push_back(entity);
  return czml;
}

int main() {
  //   Real t0_utc = GregorianToTime("2025-01-01T12:00:00");
  //   Real t0_tai = UtcToTai(t0_utc);
  //   Real tf_tai = t0_tai + 12 * SECS_HOUR;
  //   VecX tspan = VecX::LinSpaced(100, t0_tai, tf_tai);

  //   // //
  //   ************************************************************************************
  //   // // Moon Spacecraft

  //   // State coe_mop
  //   //     = ClassicalOE({6541.4, 0.6, 65.5 * RAD, 0 * RAD, 90 * RAD, 180 *
  //   RAD}, Frame::MOON_OP);
  //   // State rv_mci = ConvertFrame(t0_tai, ClassicalToCart(coe_mop, GM_MOON),
  //   Frame::MOON_CI);

  //   // auto sc_dyn = MakePtr<NBodyDynamics>();
  //   // sc_dyn->SetFrame(Frame::MOON_CI);
  //   // sc_dyn->AddBody(Body::Earth());
  //   // sc_dyn->AddBody(Body::Moon());
  //   // sc_dyn->SetTimeStep(10.0);  // [s]

  //   // GNSS spacecraft
  //   State coe = ClassicalOE({20200, 0.1, 55 * RAD, 0 * RAD, 90 * RAD, 180 *
  //   RAD}, Frame::ECI); State rv_eci = ConvertFrame(t0_tai,
  //   ClassicalToCart(coe, GM_EARTH), Frame::ECI);

  //   KeplerianDynamics<ClassicalOE> sc_dyn(GM_EARTH);
  //   MatX6 coes = sc_dyn.Propagate(coe, t0_tai, tspan);  // [s]
  //   MatX6 rvs = ClassicalToCart(coes, GM_EARTH);        // Convert to
  //   Cartesian

  //   //
  //   ********************************************************************************
  //   // Visualize with Cesium

  //   // MatX6 rvs = sc_dyn->Propagate(rv_mci, t0_tai, tspan);

  //   // Create CZML file for Cesium visualization
  //   json czml = CreateCzml(rvs, tspan);
  //   std::filesystem::path static_path = GetDataPath() / "cesium" / "static";

  //   // Save CZML to file
  //   auto czml_file = OpenFile<std::ofstream>(static_path / "satellite.czml");
  //   czml_file << czml.dump(2);
  //   czml_file.close();

  //   crow::SimpleApp app;
  //   app.loglevel(crow::LogLevel::Warning);

  //   app.route_dynamic("/")([static_path]() {
  //     std::ifstream file(static_path / "index.html");
  //     if (!file) return crow::response(404);
  //     std::ostringstream contents;
  //     contents << file.rdbuf();
  //     return crow::response(contents.str());
  //   });

  //   app.route_dynamic("/static/<string>")(
  //       [static_path](const crow::request &req, crow::response &res,
  //       std::string file) {
  //         std::ifstream in(static_path / file);
  //         if (!in) {
  //           res.code = 404;
  //           res.end();
  //           return;
  //         }
  //         std::ostringstream oss;
  //         oss << in.rdbuf();
  //         res.write(oss.str());
  //         res.end();
  //       });

  Logger::Info("CZML/server setup template is commented out in this example.");
  return 0;
}
