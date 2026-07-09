#include "lupnt/lupnt.h"
#include "py_pybind11.h"
namespace py = pybind11;

void InitDem(py::module& m) {
  m.def("load_tiff", &LoadTiff, py::arg("path"), py::arg("xlims"), py::arg("ylims"),
        py::arg("max_res"));

  // ---- LOLA / PGDA product-78 south-pole DEM loader --------------------------
  py::class_<LolaSite>(m, "LolaSite")
      .def_readonly("id", &LolaSite::id)
      .def_readonly("name", &LolaSite::name)
      .def_readonly("lat_deg", &LolaSite::lat_deg)
      .def_readonly("lon_deg", &LolaSite::lon_deg)
      .def("__repr__",
           [](const LolaSite& s) { return "<LolaSite " + s.id + " '" + s.name + "'>"; });

  py::class_<LunarDem>(m, "LunarDem")
      .def("get_elevation", &LunarDem::GetElevation, py::arg("x"), py::arg("y"),
           "Bilinearly interpolate terrain elevation [m] at native projected (x, y) [m].")
      .def("get_elevation_lat_lon", &LunarDem::GetElevationLatLon, py::arg("lat_deg"),
           py::arg("lon_deg"), "Approximate terrain elevation [m] at a latitude/longitude [deg].")
      .def_property_readonly("x", &LunarDem::x, "Pixel-center x-coordinate grid [m].")
      .def_property_readonly("y", &LunarDem::y, "Pixel-center y-coordinate grid [m].")
      .def_property_readonly("elevation", &LunarDem::elevation, "Elevation grid [m].")
      .def_property_readonly("site", &LunarDem::site)
      .def_property_readonly("rows", &LunarDem::rows)
      .def_property_readonly("cols", &LunarDem::cols)
      .def_property_readonly("x_min", &LunarDem::x_min)
      .def_property_readonly("x_max", &LunarDem::x_max)
      .def_property_readonly("y_min", &LunarDem::y_min)
      .def_property_readonly("y_max", &LunarDem::y_max)
      .def_property_readonly("center_x", &LunarDem::center_x)
      .def_property_readonly("center_y", &LunarDem::center_y);

  m.def("get_lola_sites", &GetLolaSites,
        "List the PGDA product-78 LOLA 5 m/pixel south-pole DEM sites.");
  m.def("select_lola_site", &SelectLolaSite, py::arg("lat_deg"), py::arg("lon_deg"),
        py::return_value_policy::reference,
        "Nearest PGDA product-78 site to a query latitude/longitude [deg].");
  m.def("lola_dem_url", &LolaDemUrl, py::arg("site_id"),
        "NASA PGDA download URL for a site's 5 m/pixel surface DEM GeoTIFF.");
  m.def(
      "download_lola_dem",
      [](const std::string& site_id) { return DownloadLolaDem(site_id).string(); },
      py::arg("site_id"),
      "Download (and cache under LUPNT_DATA_PATH) a site's DEM GeoTIFF; return its path.");
  m.def("load_lola_dem", &LoadLolaDem, py::arg("lat_deg"), py::arg("lon_deg"),
        py::arg("half_width_m") = 5000.0, py::arg("max_res") = 20.0,
        "Load the appropriate NASA PGDA DEM for a query latitude/longitude [deg].");
}
