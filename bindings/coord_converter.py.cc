#include <lupnt/physics/coord_converter.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_coord_converter(py::module &m) {
  py::enum_<CoordSystem>(m, "CoordSystem")
      .value("ITRF", CoordSystem::ITRF)
      .value("GCRF", CoordSystem::GCRF)
      .value("ICRF", CoordSystem::ICRF)
      .value("SER", CoordSystem::SER)
      .value("GSE", CoordSystem::GSE)
      .value("EME", CoordSystem::EME)
      .value("MOD", CoordSystem::MOD)
      .value("TOD", CoordSystem::TOD)
      .value("EMR", CoordSystem::EMR)
      .value("MI", CoordSystem::MI)
      .value("PA", CoordSystem::PA)
      .value("ME", CoordSystem::ME)
      .value("RTN", CoordSystem::RTN)
      .value("CoordSystemCount", CoordSystem::CoordSystemCount)
      .value("NONE", CoordSystem::NONE)
      .export_values();

  py::class_<CoordConverter>(m, "CoordConverter")
      .def_static(
          "convert",
          [](const VectorX &rv_in, const real epoch,
             const CoordSystem coord_sys_in, const CoordSystem coord_sys_out) {
            return CoordConverter::Convert(rv_in, epoch, coord_sys_in,
                                           coord_sys_out);
          },
          "Convert frame (Frame ID input)");
}