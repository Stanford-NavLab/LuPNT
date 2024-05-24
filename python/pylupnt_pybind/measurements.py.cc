#include <lupnt/lupnt.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace lupnt;

void init_measurements(py::module &m) {
  m.def(
      "compute_occultation",
      [](double epoch, const Vector3d &r1, const Vector3d &r2, Frame cs1,
         Frame cs2, const std::vector<NaifId> &bodies) -> VectorXd {
        return Occultation::ComputeOccultation(epoch, r1.cast<real>().eval(),
                                               r2.cast<real>().eval(), cs1, cs2,
                                               bodies)
            .cast<double>();
      },
      "Compute occultation between two points", py::arg("epoch"), py::arg("r1"),
      py::arg("r2"), py::arg("cs1"), py::arg("cs2"), py::arg("bodies"));
  m.def(
      "compute_occultation",
      [](double epoch, Matrixd<-1, 3> &r1, Matrixd<-1, 3> &r2, Frame cs1,
         Frame cs2, const std::vector<NaifId> &bodies) -> MatrixXd {
        return Occultation::ComputeOccultation(epoch, r1.cast<real>().eval(),
                                               r2.cast<real>().eval(), cs1, cs2,
                                               bodies)
            .cast<double>();
      },
      "Compute occultation between two points", py::arg("epoch"), py::arg("r1"),
      py::arg("r2"), py::arg("cs1"), py::arg("cs2"), py::arg("bodies"));
  m.def(
      "compute_occultation",
      [](const VectorXd &epoch, const Matrixd<-1, 3> &r1,
         const Matrixd<-1, 3> &r2, Frame cs1, Frame cs2,
         const std::vector<NaifId> &bodies) -> MatrixXd {
        return Occultation::ComputeOccultation(
                   epoch.cast<real>().eval(), r1.cast<real>().eval(),
                   r2.cast<real>().eval(), cs1, cs2, bodies)
            .cast<double>();
      },
      "Compute occultation between two points", py::arg("epoch"), py::arg("r1"),
      py::arg("r2"), py::arg("cs1"), py::arg("cs2"), py::arg("bodies"));
}
