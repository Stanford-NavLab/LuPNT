#include <lupnt/interfaces/eop.h>

#include <filesystem>
#include <string>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitEop(py::module& m) {
  m.def(
      "load_eop_file_data",
      [](const std::string& filepath, bool force) {
        LoadEopFileData(std::filesystem::path(filepath), force);
      },
      py::arg("filepath"), py::arg("force") = false,
      "Parse an IERS EOP 14 C04 (IAU1980) file at `filepath` into the global EOP table. "
      "No-op if already loaded unless `force`.");

  m.def("load_latest_eop_from_iers", &LoadLatestEopFromIers, py::arg("force") = true,
        "Download the latest IERS EOP 14 C04 (IAU1980) series and load it; returns False (falling "
        "back to bundled data) on download failure.");

  m.def(
      "get_eop_data",
      [](Real mjd_utc) {
        EopData data = GetEopData(mjd_utc);
        py::dict d;
        d["x_pole"] = data.x_pole.val();
        d["y_pole"] = data.y_pole.val();
        d["ut1_utc"] = data.ut1_utc.val();
        d["lod"] = data.lod.val();
        d["dpsi"] = data.dpsi.val();
        d["deps"] = data.deps.val();
        d["dx_pole"] = data.dx_pole.val();
        d["dy_pole"] = data.dy_pole.val();
        d["tai_utc"] = data.tai_utc.val();
        return d;
      },
      py::arg("mjd_utc"),
      "Interpolate Earth orientation parameters at `mjd_utc` [MJD, UTC] from the loaded IERS EOP "
      "table, as a dict (x_pole/y_pole [rad], ut1_utc [s], lod [s], dpsi/deps [rad], "
      "dx_pole/dy_pole [rad], tai_utc [s]).");

  m.def(
      "get_eop_file_data",
      []() {
        EopFileData* file_data = GetEopFileData();
        py::dict d;
        d["mjds_utc"] = file_data->mjds_utc;
        d["x"] = file_data->x;
        d["y"] = file_data->y;
        d["ut1_utc"] = file_data->ut1_utc;
        d["lod"] = file_data->lod;
        d["dpsi"] = file_data->dpsi;
        d["deps"] = file_data->deps;
        d["x_err"] = file_data->xErr;
        d["y_err"] = file_data->yErr;
        d["ut1_utc_err"] = file_data->ut1_utc_err;
        d["lod_err"] = file_data->lod_err;
        d["dpsi_err"] = file_data->dpsi_err;
        d["deps_err"] = file_data->deps_err;
        return d;
      },
      "Return the full raw loaded IERS EOP table (time series + formal errors) as a dict of "
      "arrays, loading the bundled file first if none is loaded yet.");
}
