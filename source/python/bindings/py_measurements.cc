#include <lupnt/lupnt.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace lupnt;

void init_measurements(py::module &m) {
  m.def(
      "compute_occultation",
      [](double epoch, const Vec3d &r1, const Vec3d &r2, Frame cs1, Frame cs2,
         const std::vector<NaifId> &bodies) -> VecXd {
        return Occultation::ComputeOccultation(epoch, r1.cast<Real>().eval(),
                                               r2.cast<Real>().eval(), cs1, cs2, bodies)
            .cast<double>();
      },
      "Compute occultation between two points", py::arg("epoch"), py::arg("r1"), py::arg("r2"),
      py::arg("cs1"), py::arg("cs2"), py::arg("bodies"));
  m.def(
      "compute_occultation",
      [](double epoch, Matd<-1, 3> &r1, Matd<-1, 3> &r2, Frame cs1, Frame cs2,
         const std::vector<NaifId> &bodies) -> VecXd {
        return Occultation::ComputeOccultation(epoch, r1.cast<Real>().eval(),
                                               r2.cast<Real>().eval(), cs1, cs2, bodies)
            .cast<double>();
      },
      "Compute occultation between two points", py::arg("epoch"), py::arg("r1"), py::arg("r2"),
      py::arg("cs1"), py::arg("cs2"), py::arg("bodies"));
  m.def(
      "compute_occultation",
      [](const VecXd &epoch, const Matd<-1, 3> &r1, const Matd<-1, 3> &r2, Frame cs1, Frame cs2,
         const std::vector<NaifId> &bodies) -> VecXd {
        return Occultation::ComputeOccultation(epoch.cast<Real>().eval(), r1.cast<Real>().eval(),
                                               r2.cast<Real>().eval(), cs1, cs2, bodies)
            .cast<double>();
      },
      "Compute occultation between two points", py::arg("epoch"), py::arg("r1"), py::arg("r2"),
      py::arg("cs1"), py::arg("cs2"), py::arg("bodies"));
}
