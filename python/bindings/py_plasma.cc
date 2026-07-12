/**
 * @file py_plasma.cc
 * @brief Python bindings for the plasma (pecsim) module — ionospheric/
 *        plasmaspheric electron density and GNSS ray-tracing.
 *
 * Adapted from pecsim/src/python/bindings/ for integration into LuPNT.
 * The C++ implementation lives in cpp/lupnt/environment/plasma/ and uses
 * namespace pecsim::.  Data files are loaded from PECSIMPY_BASE_PATH
 * (set to $LUPNT_DATA_PATH/plasma by the pixi activation env).
 */

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "lupnt/environment/plasma/plasma.h"

namespace py = pybind11;
using namespace pecsim;

// ---- bind_constants -------------------------------------------------------

static void bind_constants(py::module& m) {
  m.attr("PLASMA_BASE_PATH") = get_base_path();
  m.attr("RE_PLASMA") = py::float_(RE);
  m.attr("C_PLASMA") = py::float_(C);
  m.attr("PI_PLASMA") = py::float_(PI);
  m.attr("SECS_DAY_PLASMA") = py::float_(SECS_DAY);
  m.attr("RAD2DEG_PLASMA") = py::float_(RAD2DEG);
  m.attr("DEG2RAD_PLASMA") = py::float_(DEG2RAD);
  m.attr("TECU_PLASMA") = py::float_(TECU);
  m.attr("GM_EARTH_PLASMA") = py::float_(GM_EARTH);
  m.attr("freq_L1") = py::float_(freq_L1);
  m.attr("freq_L2") = py::float_(freq_L2);
  m.attr("freq_L5") = py::float_(freq_L5);
}

// ---- bind_core ------------------------------------------------------------

static void bind_core(py::module& m) {
  m.def(
      "get_plasma_base_path", []() { return get_base_path(); },
      "Get the base path used to locate plasma runtime data (iri/, kp/).");
  m.def(
      "set_plasma_base_path", [](const std::string& path) { set_base_path(path); }, py::arg("path"),
      "Override the base path used to locate plasma runtime data.");
}

// ---- bind_time ------------------------------------------------------------

static void bind_time(py::module& m) {
  py::class_<DateTime>(m, "DateTime",
                       "Calendar date/time as year, day-of-year, hour, minute, second.")
      .def(py::init<int, int, int, int, double>(), py::arg("year"), py::arg("doy"), py::arg("hour"),
           py::arg("min"), py::arg("sec"))
      .def_readwrite("year", &DateTime::year, "Year.")
      .def_readwrite("doy", &DateTime::doy, "Day of year.")
      .def_readwrite("hour", &DateTime::hour, "Hour [0-23].")
      .def_readwrite("min", &DateTime::min, "Minute [0-59].")
      .def_readwrite("sec", &DateTime::sec, "Second [0-60).");

  m.def("datetime_to_itime", &datetime_to_itime, py::arg("datetime"),
        "Convert a DateTime to an itime array [year+doy, milliseconds-of-day].");
  m.def("itime_to_datetime", &itime_to_datetime, py::arg("itime"),
        "Convert an itime array [year+doy, milliseconds-of-day] to a DateTime.");
  m.def("mjd_to_datetime", &mjd_to_datetime, py::arg("mjd"),
        "Convert a Modified Julian Date to a DateTime.");
  m.def("datetime_to_mjd", &datetime_to_mjd, py::arg("datetime"),
        "Convert a DateTime to a Modified Julian Date.");
  m.def("gregorian_to_mjd", &gregorian_to_mjd, py::arg("year"), py::arg("month"), py::arg("day"),
        py::arg("hour") = 0, py::arg("min") = 0, py::arg("sec") = 0.0,
        "Convert a Gregorian calendar date to a Modified Julian Date.");
  m.def("tj2000_to_mjd", &tj2000_to_mjd, py::arg("t_j2000"),
        "Convert seconds past J2000 to a Modified Julian Date.");
  m.def("mjd_to_tj2000", &mjd_to_tj2000, py::arg("mjd"),
        "Convert a Modified Julian Date to seconds past J2000.");
  m.def("long_to_lt", &long_to_lt, py::arg("longitude"),
        "Convert longitude [rad] to local time [hours].");
  m.def("lt_to_long", &lt_to_long, py::arg("amlt"),
        "Convert magnetic local time [hours] to longitude [rad].");
}

// ---- bind_gcpm ------------------------------------------------------------

static void set_iri_model_(const std::string& model_name) {
  if (model_name == "IRI2007") {
    set_iri_model(IRIModel::IRI_2007);
  } else if (model_name == "IRI2020") {
    set_iri_model(IRIModel::IRI_2020);
  } else {
    throw std::invalid_argument("Invalid IRI model name. Available: IRI2007, IRI2020.");
  }
}

static void bind_gcpm(py::module& m) {
  py::class_<IRI2007Option>(m, "IRI2007Option", "Configuration flags for the IRI-2007 model.")
      .def(py::init<>())
      .def_readwrite("R12", &IRI2007Option::R12,
                     "R12 sunspot index; >0 uses the value, -1 historical/projected (with storm "
                     "model), -2 without storm model.")
      .def("update_jf_2007", &IRI2007Option::update_jf_2007,
           "Rebuild the internal IRI-2007 jf option flags from the current fields.");

  py::class_<IRI2020Option>(m, "IRI2020Option", "Configuration flags for the IRI-2020 model.")
      .def(py::init<>())
      .def_readwrite("compute_teti", &IRI2020Option::compute_teti,
                     "Compute electron/ion temperatures.")
      .def_readwrite("compute_ni", &IRI2020Option::compute_ni, "Compute ion densities.")
      .def_readwrite("output_messages", &IRI2020Option::output_messages,
                     "Emit IRI status messages.")
      .def_readwrite("output_to_text", &IRI2020Option::output_to_text, "Write IRI output to text.")
      .def_readwrite("plasma_model", &IRI2020Option::plasma_model,
                     "Plasmasphere model (1: Ozhogin, 0: Gallagher).")
      .def_readwrite("without_plasmapause", &IRI2020Option::without_plasmapause,
                     "Disable the plasmapause when set.")
      .def_readwrite("R12", &IRI2020Option::R12, "R12 sunspot index (unused by IRI-2020).")
      .def("update_jf_2020", &IRI2020Option::update_jf_2020,
           "Rebuild the internal IRI-2020 jf option flags from the current fields.");

  m.def("set_iri_model", &set_iri_model_, py::arg("model_name") = "IRI2007",
        "Select the IRI ionosphere model (\"IRI2007\" or \"IRI2020\").");
  m.def("get_iri_model", &get_iri_model_str, "Name of the currently selected IRI model.");
  m.def("set_iri2007_option", &set_iri2007_option, py::arg("option"),
        "Set the global IRI-2007 configuration.");
  m.def("set_iri2020_option", &set_iri2020_option, py::arg("option"),
        "Set the global IRI-2020 configuration.");
  m.def("get_kp_index", &get_kp_index, py::arg("datetime"),
        "Kp geomagnetic activity index at a given DateTime.");
  m.def("gcpm_v24", &gcpm_v24, py::arg("datetime"), py::arg("r_RE"), py::arg("amlt"),
        py::arg("alatr"), py::arg("akp") = -1.0, "C++ GCPM v2.4: returns [ne_cm3, H+, He+, O+].");
  m.def("gcpm_v24_fortran", &gcpm_v24_fortran, py::arg("datetime"), py::arg("r_RE"),
        py::arg("amlt"), py::arg("alatr"), py::arg("akp") = -1.0,
        "Fortran GCPM v2.4: returns [ne_cm3, H+, He+, O+].");
}

// ---- bind_orbit -----------------------------------------------------------

static void bind_orbit(py::module& m) {
  py::class_<Satellite>(m, "Satellite",
                        "Keplerian satellite propagated with Earth two-body dynamics (km units).")
      .def(py::init<int, const Vec6d&, double, double>(), py::arg("id"), py::arg("coe"),
           py::arg("epoch_utc"), py::arg("GM") = GM_EARTH)
      .def("propagate", &Satellite::propagate, py::arg("epoch_utc_new"),
           "Propagate the satellite state to a new UTC epoch [s past J2000].")
      .def("get_pos", &Satellite::get_pos, "Position in ECI coordinates [km].")
      .def("get_pos_epoch", &Satellite::get_pos_epoch, py::arg("dt"),
           "Position [km] propagated forward by dt [s] from the current epoch.")
      .def_readwrite("id_", &Satellite::id_, "Satellite ID (PRN).")
      .def_readwrite("posvel_", &Satellite::posvel_, "Position/velocity state [km, km/s].")
      .def_readwrite("coe_", &Satellite::coe_, "Classical orbital elements (a, e, i, Omega, w, M).")
      .def_readwrite("epoch_utc_", &Satellite::epoch_utc_, "Epoch [s UTC past J2000].")
      .def_readwrite("GM_", &Satellite::GM_, "Gravitational parameter [km^3/s^2].");

  m.def("wrap2pi", &wrap2pi, py::arg("angle"), "Wrap an angle [rad] to [0, 2*pi).");
  m.def("coe2cart", &coe2cart, py::arg("coe"), py::arg("GM"),
        "Convert classical orbital elements to a Cartesian state given GM [km^3/s^2].");
  m.def("cart2coe", &cart2coe, py::arg("posvel"), py::arg("GM"),
        "Convert a Cartesian state to classical orbital elements given GM [km^3/s^2].");
  m.def("propagate_coe", &propagate_coe, py::arg("coe"), py::arg("GM"), py::arg("dt"),
        "Propagate classical orbital elements by dt [s] under two-body motion.");
  m.def("compute_vis", &compute_vis, py::arg("pos1"), py::arg("pos2"), py::arg("radius"),
        py::arg("max_angle"),
        "True if pos1 and pos2 [km] have line-of-sight past a blocking sphere of given radius [km] "
        "and within max_angle [rad].");
  m.def("compute_min_altitude", &compute_min_altitude, py::arg("pos1"), py::arg("pos2"),
        py::arg("radius"),
        "Minimum altitude [km] of the line segment pos1-pos2 above a sphere of given radius [km].");
  m.def("mean2true", &mean2true, py::arg("M"), py::arg("e"),
        "Mean anomaly [rad] to true anomaly [rad] for eccentricity e.");
  m.def("true2mean", &true2mean, py::arg("nu"), py::arg("e"),
        "True anomaly [rad] to mean anomaly [rad] for eccentricity e.");
  m.def("ecc2true", &ecc2true, py::arg("E"), py::arg("e"),
        "Eccentric anomaly [rad] to true anomaly [rad] for eccentricity e.");
  m.def("mean2ecc", &mean2ecc, py::arg("M"), py::arg("e"),
        "Mean anomaly [rad] to eccentric anomaly [rad] for eccentricity e.");
  m.def("solve_lt", &solve_lt, py::arg("sat"), py::arg("rx_pos"), py::arg("epoch_utc_rx"),
        "Solve light time to find the transmitter position [km] seen by a receiver at rx_pos.");
  m.def("setup_gnss_constellation", &setup_gnss_constellation, py::arg("filename"),
        "Build a list of Satellites for a GNSS constellation from a definition file.");
}

// ---- bind_tec -------------------------------------------------------------

static void bind_tec(py::module& m) {
  py::class_<RayTraceConfig>(m, "RayTraceConfig", "Configuration for ionospheric ray tracing.")
      .def(py::init<>())
      .def_readwrite("freq_Hz", &RayTraceConfig::freq_Hz, "Signal frequency [Hz].")
      .def_readwrite("step_size", &RayTraceConfig::step_size, "Ray-tracing step size [km].")
      .def_readwrite("correction", &RayTraceConfig::correction, "Apply endpoint correction.")
      .def_readwrite("fine_correction", &RayTraceConfig::fine_correction, "Apply fine correction.")
      .def_readwrite("cutoff_r", &RayTraceConfig::cutoff_r, "Outer cutoff radius for tracing [km].")
      .def_readwrite("gradn_dx", &RayTraceConfig::gradn_dx,
                     "Finite-difference step for the refractive-index gradient [km].")
      .def_readwrite("integ_method", &RayTraceConfig::integ_method,
                     "Integration method (\"RK4\" or \"Euler\").")
      .def_readwrite("correction_method", &RayTraceConfig::correction_method,
                     "Correction method (\"grid\" or \"newton\").")
      .def_readwrite("kp", &RayTraceConfig::kp, "Kp index for the ionosphere model.")
      .def_readwrite("rz12", &RayTraceConfig::rz12,
                     "IRI R12 sunspot index; >0 uses value, -1 with storm model, -2 without.")
      .def_readwrite("use_fortran_gcpm", &RayTraceConfig::use_fortran_gcpm,
                     "Use the Fortran GCPM implementation.")
      .def_readwrite("corr_tol", &RayTraceConfig::corr_tol, "Endpoint correction tolerance [m].")
      .def_readwrite("compute_higher_order", &RayTraceConfig::compute_higher_order,
                     "Compute second-order (and higher) delays.")
      .def_readwrite("use_adaptive_step", &RayTraceConfig::use_adaptive_step,
                     "Use adaptive step-size control.")
      .def_readwrite("straight_ray", &RayTraceConfig::straight_ray, "Assume a straight ray path.");

  py::class_<PathProfile>(m, "PathProfile",
                          "Result of tracing a ray: per-step samples plus summary quantities.")
      .def(py::init<>())
      .def_readwrite("s", &PathProfile::s, "Per-step path length [km].")
      .def_readwrite("tec_section", &PathProfile::tec_section,
                     "Per-step total electron content [TECU].")
      .def_readwrite("r", &PathProfile::r, "Per-step geocentric radial distance [km].")
      .def_readwrite("pos_eci", &PathProfile::pos_eci, "Per-step position in ECI [km].")
      .def_readwrite("az_dir", &PathProfile::az_dir, "Per-step azimuth direction [rad].")
      .def_readwrite("el_dir", &PathProfile::el_dir, "Per-step elevation direction [rad].")
      .def_readwrite("dist_to_line", &PathProfile::dist_to_line,
                     "Per-step distance to the straight line-of-sight [km].")
      .def_readwrite("sf", &PathProfile::sf, "Final ray distance [km].")
      .def_readwrite("dir_start", &PathProfile::dir_start, "Direction vector at the start.")
      .def_readwrite("dir_end", &PathProfile::dir_end, "Direction vector at the end.")
      .def_readwrite("final_pos", &PathProfile::final_pos,
                     "Final position in Cartesian coordinates [km].")
      .def_readwrite("corr_final_pos_err", &PathProfile::corr_final_pos_err,
                     "Final position error [km] from insufficient correction.")
      .def_readwrite("corr_final_time_err", &PathProfile::corr_final_time_err,
                     "Final time error [s] from insufficient correction.")
      .def_readwrite("t_tx", &PathProfile::t_tx, "Transmission epoch [s past J2000].")
      .def_readwrite("t_rx", &PathProfile::t_rx, "Reception epoch [s past J2000].")
      .def_readwrite("prop_time_total", &PathProfile::prop_time_total,
                     "Total propagation time [s].")
      .def_readwrite("dist_bend_m", &PathProfile::dist_bend_m,
                     "Extra path length from bending [m].")
      .def_readwrite("dist_straight_km", &PathProfile::dist_straight_km,
                     "Straight-line distance [km].")
      .def_readwrite("total_delay_m", &PathProfile::total_delay_m, "Total group delay [m].")
      .def_readwrite("tec_delay_m", &PathProfile::tec_delay_m, "First-order TEC delay [m].")
      .def_readwrite("tec_delay_bend_m", &PathProfile::tec_delay_bend_m,
                     "Additional TEC delay from ray bending [m].")
      .def_readwrite("second_delay_m", &PathProfile::second_delay_m, "Second-order delay [m].")
      .def_readwrite("third_delay_m", &PathProfile::third_delay_m, "Third-order delay [m].")
      .def_readwrite("tecu", &PathProfile::tecu, "Total electron content [TECU].")
      .def_readwrite("tecu_bend", &PathProfile::tecu_bend, "TEC contribution from bending [TECU].")
      .def_readwrite("max_sep_line_m", &PathProfile::max_sep_line_m,
                     "Maximum separation from the straight line [m].");

  m.def(
      "trace_ray",
      [](double epoch, const Eigen::Vector3d& pos_tx, const Eigen::Vector3d& pos_rx,
         const RayTraceConfig& config, bool debug_prop, bool debug_corr) {
        py::scoped_ostream_redirect stream(std::cout, py::module_::import("sys").attr("stdout"));
        return trace_ray(epoch, pos_tx, pos_rx, config, debug_prop, debug_corr);
      },
      py::arg("epoch"), py::arg("pos_tx"), py::arg("pos_rx"), py::arg("config"),
      py::arg("debug_prop") = false, py::arg("debug_corr") = false,
      "Trace a ray from pos_tx to pos_rx and return a PathProfile.");

  m.def("compute_ne", &compute_ne, py::arg("t_j2000"), py::arg("pos_geo"), py::arg("config"),
        py::arg("debug") = false, "Electron density [cm^-3] at a given position and time.");
  m.def("get_iono_params", &get_iono_params, py::arg("t_j2000"), py::arg("kp") = -1.0,
        "Returns [f107, rz12, hmf2_km] ionospheric parameters.");
  m.def("compute_B", &compute_B, py::arg("t_j2000"), py::arg("pos_geo"), py::arg("config"),
        py::arg("debug") = false, "Magnetic field vector [nT] at a given ECEF position.");
  m.def("refractive_index_neB", &refractive_index_neB, py::arg("ne_m3"), py::arg("freq_Hz"),
        py::arg("B") = 0.0, py::arg("cos_theta") = 0.0, py::arg("compute_higher_order") = true,
        "Refractive index of the ionosphere for given ne and frequency.");
}

// ---- InitPlasma -----------------------------------------------------------

void InitPlasma(py::module& m) {
  py::add_ostream_redirect(m, "ostream_redirect_plasma");
  bind_constants(m);
  bind_core(m);
  bind_time(m);
  bind_gcpm(m);
  bind_orbit(m);
  bind_tec(m);
}
