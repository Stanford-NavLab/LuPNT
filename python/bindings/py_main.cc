#define PYBIND11_DETAILED_ERROR_MESSAGES
#include <pybind11/pybind11.h>

namespace py = pybind11;

void InitAutodiff(py::module& m);
void InitConstants(py::module& m);
void InitFrameConverter(py::module& m);
void InitTimeConverter(py::module& m);
void InitDynamics(py::module& m);
void InitScUdFilter(py::module& m);
void InitConversions(py::module& m);
void InitAgents(py::module& m);
void InitAntenna(py::module& m);
void InitTle(py::module& m);
void InitKernels(py::module& m);
void InitMathUtils(py::module& m);
void InitFile(py::module& m);
void InitConfig(py::module& m);
void InitLogger(py::module& m);
void InitDem(py::module& m);
void InitPlasma(py::module& m);
void InitGnss(py::module& m);
void InitEop(py::module& m);
void InitIslOdts(py::module& m);
void InitSimulation(py::module& m);
void InitGroundStationOdts(py::module& m);
void InitEphemeris(py::module& m);
void InitGnssOdts(py::module& m);
void InitSurfaceNav(py::module& m);
void InitLanderNav(py::module& m);
void InitEpoch(py::module& m);
// Defined in py_spice_interface.cc (snake_case to match that file).
void init_spice_interface(py::module& m);
// Defined in py_forces.cc (snake_case to match that file). Exposes the
// free-function acceleration models (gravity field, point mass, SRP,
// relativistic correction, shadow/illumination) for scripting and for
// cross-validating force models from Python.
void init_forces(py::module& m);

PYBIND11_MODULE(_pylupnt, m) {
  InitAutodiff(m);
  InitConstants(m);
  InitTimeConverter(m);
  InitFrameConverter(m);
  InitDynamics(m);
  InitScUdFilter(m);
  InitConversions(m);
  // InitAgents(m);
  InitAntenna(m);
  InitTle(m);
  InitKernels(m);
  InitMathUtils(m);
  InitFile(m);
  InitConfig(m);
  InitLogger(m);
  InitDem(m);
  InitPlasma(m);
  InitGnss(m);
  InitEop(m);
  InitSimulation(m);
  InitIslOdts(m);
  InitGroundStationOdts(m);
  InitEphemeris(m);
  InitGnssOdts(m);
  InitSurfaceNav(m);
  InitLanderNav(m);
  InitEpoch(m);
  // Exposes the `pnt.spice` submodule (kernel loading, SPICE-based frame/time
  // conversions, and the scalar time-ephemeris accessor used to read the
  // DE440t TT-TDB and the lunar TDB-TCL / TDB-TL kernels).
  init_spice_interface(m);
  init_forces(m);
}
