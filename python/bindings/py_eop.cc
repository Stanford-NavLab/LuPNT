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
        "back to bundled data) on download failure. Note C04 is a retrospective final series that "
        "lags by months -- for present-day or future epochs use load_latest_eop_finals_from_iers "
        "instead.");

  py::enum_<EopSource>(m, "EopSource", "Which IERS product the loaded EOP table came from.")
      .value("C04", EopSource::C04, "IERS EOP 14 C04 final series; retrospective, no predictions.")
      .value("Finals", EopSource::Finals,
             "IERS finals.all (Bulletin A); rapid solution plus ~1 year of predictions.");

  m.def(
      "set_eop_source",
      [](EopSource source, const std::string& filepath) {
        SetEopSource(source, std::filesystem::path(filepath));
      },
      py::arg("source"), py::arg("filepath") = "",
      "Declare which EOP product to use, resolved when the table is first needed. Call order does "
      "not matter: if nothing is loaded yet the choice is applied on first use, and if a "
      "different table is already loaded it is replaced immediately -- unlike a bare "
      "load_*(force=False), which is silently ignored once anything has touched Earth "
      "orientation. Defaults to C04 (reproducible but retrospective); select Finals for any run "
      "at a present-day or future epoch. Raises if the file does not exist.");

  m.def("get_eop_source", &GetEopSource,
        "The currently declared EOP source, which may not be loaded yet. For what is actually "
        "loaded, use get_eop_coverage()['source'].");

  py::enum_<EopNutation>(m, "EopNutation",
                         "Meaning of a file's celestial-pole columns. The finals variants share a "
                         "layout but not a meaning, so this cannot be auto-detected.")
      .value("Iau1980", EopNutation::Iau1980, "dPsi, dEps (C04, finals.all).")
      .value("Iau2000", EopNutation::Iau2000,
             "dX, dY -- CIP offsets to IAU 2006/2000A (finals2000A.all).");

  py::class_<EopPerturbation>(
      m, "EopPerturbation",
      "Additive EOP perturbation (perturbed = nominal + delta), injected inside get_eop_data so "
      "it reaches polar motion, UT1/sidereal rotation, LOD and the CIP offsets at once.")
      .def(py::init<>())
      .def_readwrite("dx_pole", &EopPerturbation::dx_pole, "Polar motion x [rad]")
      .def_readwrite("dy_pole", &EopPerturbation::dy_pole, "Polar motion y [rad]")
      .def_readwrite("dut1", &EopPerturbation::dut1, "UT1-UTC [s]")
      .def_readwrite("dlod", &EopPerturbation::dlod, "Length of day [s]")
      .def_readwrite("ddX", &EopPerturbation::ddX, "CIP offset dX [rad]")
      .def_readwrite("ddY", &EopPerturbation::ddY, "CIP offset dY [rad]");

  m.def("set_eop_perturbation", &SetEopPerturbation, py::arg("perturbation"),
        "Apply a constant EOP perturbation to every subsequent EOP evaluation.");
  m.def("set_eop_perturbation_function", &SetEopPerturbationFunction, py::arg("fn"),
        "Apply a time-varying EOP perturbation; `fn(mjd_utc)` must be thread-safe.");
  m.def("clear_eop_perturbation", &ClearEopPerturbation, "Remove any EOP perturbation.");
  m.def("get_eop_perturbation", &GetEopPerturbation, py::arg("mjd_utc"),
        "The perturbation in force at `mjd_utc` (all-zero when none is set).");
  m.def("eop_has_celestial_pole_offsets", &EopHasCelestialPoleOffsets,
        "True when an IAU2000 table is loaded or a perturbation is active, i.e. when the CIP "
        "offsets can be non-zero.");

  m.def(
      "load_eop_finals_file_data",
      [](const std::string& filepath, bool force, EopNutation nutation) {
        LoadEopFinalsFileData(std::filesystem::path(filepath), force, nutation);
      },
      py::arg("filepath"), py::arg("force") = false, py::arg("nutation") = EopNutation::Iau1980,
      "Parse an IERS finals ('Bulletin A', IAU1980) file at `filepath` into the global EOP table. "
      "Unlike C04 this reaches the present day and carries ~1 year of predictions. No-op if data "
      "is already loaded unless `force`.");

  m.def("load_latest_eop_finals_from_iers", &LoadLatestEopFinalsFromIers, py::arg("force") = true,
        "Download the latest IERS finals.all (Bulletin A, IAU1980) file and load it; returns "
        "False on download failure, leaving any already-loaded table in place. This is the only "
        "loader that reaches the present day and beyond -- prefer it for any run at a present-day "
        "or future epoch.");

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
        d["sigma_x_pole"] = data.sigma_x_pole.val();
        d["sigma_y_pole"] = data.sigma_y_pole.val();
        d["sigma_ut1_utc"] = data.sigma_ut1_utc.val();
        d["dX"] = data.dX.val();
        d["dY"] = data.dY.val();
        return d;
      },
      py::arg("mjd_utc"),
      "Interpolate Earth orientation parameters at `mjd_utc` [MJD, UTC] from the loaded IERS EOP "
      "table, as a dict (x_pole/y_pole [rad], ut1_utc [s], lod [s], dpsi/deps [rad], and the C04 "
      "formal errors sigma_x_pole/sigma_y_pole [rad], sigma_ut1_utc [s]). Outside the table range "
      "the nearest endpoint is held constant -- see get_eop_coverage().");

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
        d["is_prediction"] = file_data->is_prediction;
        d["source"] = file_data->source;
        d["nutation"] = file_data->nutation;
        d["dX"] = file_data->dX;
        d["dY"] = file_data->dY;
        d["dX_err"] = file_data->dX_err;
        d["dY_err"] = file_data->dY_err;
        return d;
      },
      "Return the full raw loaded IERS EOP table (time series + formal errors) as a dict of "
      "arrays, loading the bundled file first if none is loaded yet.");

  m.def(
      "get_eop_coverage",
      []() {
        EopCoverage c = GetEopCoverage();
        py::dict d;
        d["mjd_first"] = c.mjd_first;
        d["mjd_last"] = c.mjd_last;
        d["mjd_last_measured"] = c.mjd_last_measured;
        d["source"] = c.source;
        return d;
      },
      "Return the epoch coverage of the loaded EOP table as a dict {mjd_first, mjd_last, "
      "mjd_last_measured, source} [MJD, UTC], loading the bundled file first if none is loaded "
      "yet. Past mjd_last_measured the values are predicted; past mjd_last get_eop_data() holds "
      "the nearest endpoint constant instead of extrapolating.");
}
