// lupnt
#include <lupnt/core/constants.h>
#include <lupnt/interfaces/kernels.h>
#include <lupnt/numerics/math_utils.h>

// pybind11
#include <string>

#include "py_pybind11.h"
#include "py_vectorized_macros.h"

namespace py = pybind11;
using namespace lupnt;

void InitKernels(py::module& m) {
  // GetBodyPosVel
  m.def("get_body_pos_vel", py::overload_cast<Real, BodyId, BodyId, Frame>(&GetBodyPosVel),
        py::arg("t_tdb"), py::arg("center"), py::arg("target"), py::arg("frame"));

  m.def("get_body_pos_vel",
        py::overload_cast<Real, BodyId, BodyId, Frame, const UnitSystem&, CoordinateScale>(
            &GetBodyPosVel),
        py::arg("t_tdb"), py::arg("center"), py::arg("target"), py::arg("frame"),
        py::arg("units") = SI_UNITS, py::arg("coordinate_scale") = CoordinateScale::TDB);

  m.def("get_body_pos_vel", py::overload_cast<const VecX&, BodyId, BodyId, Frame>(&GetBodyPosVel),
        py::arg("t_tdb"), py::arg("center"), py::arg("target"), py::arg("frame"));

  m.def("get_body_pos_vel",
        py::overload_cast<const VecX&, BodyId, BodyId, Frame, const UnitSystem&, CoordinateScale>(
            &GetBodyPosVel),
        py::arg("t_tdb"), py::arg("center"), py::arg("target"), py::arg("frame"),
        py::arg("units") = SI_UNITS, py::arg("coordinate_scale") = CoordinateScale::TDB);

  m.def("get_body_pos_vel", py::overload_cast<Real, BodyId, Frame>(&GetBodyPosVel),
        py::arg("t_tdb"), py::arg("target"), py::arg("frame"));

  m.def("get_body_pos_vel",
        py::overload_cast<Real, BodyId, Frame, const UnitSystem&, CoordinateScale>(&GetBodyPosVel),
        py::arg("t_tdb"), py::arg("target"), py::arg("frame"), py::arg("units") = SI_UNITS,
        py::arg("coordinate_scale") = CoordinateScale::TDB);

  m.def("get_body_pos_vel", py::overload_cast<const VecX&, BodyId, Frame>(&GetBodyPosVel),
        py::arg("t_tdb"), py::arg("target"), py::arg("frame"));

  m.def("get_body_pos_vel",
        py::overload_cast<const VecX&, BodyId, Frame, const UnitSystem&, CoordinateScale>(
            &GetBodyPosVel),
        py::arg("t_tdb"), py::arg("target"), py::arg("frame"), py::arg("units") = SI_UNITS,
        py::arg("coordinate_scale") = CoordinateScale::TDB);

  m.def("get_body_pos", py::overload_cast<Real, BodyId, BodyId, Frame>(&GetBodyPos),
        py::arg("t_tdb"), py::arg("center"), py::arg("target"), py::arg("frame"));

  m.def("get_body_pos",
        py::overload_cast<Real, BodyId, BodyId, Frame, const UnitSystem&, CoordinateScale>(
            &GetBodyPos),
        py::arg("t_tdb"), py::arg("center"), py::arg("target"), py::arg("frame"),
        py::arg("units") = SI_UNITS, py::arg("coordinate_scale") = CoordinateScale::TDB);

  m.def("get_body_pos", py::overload_cast<Real, BodyId, Frame>(&GetBodyPos), py::arg("t_tdb"),
        py::arg("target"), py::arg("frame"));

  m.def("get_body_pos",
        py::overload_cast<Real, BodyId, Frame, const UnitSystem&, CoordinateScale>(&GetBodyPos),
        py::arg("t_tdb"), py::arg("target"), py::arg("frame"), py::arg("units") = SI_UNITS,
        py::arg("coordinate_scale") = CoordinateScale::TDB);
}
