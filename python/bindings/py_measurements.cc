#include <lupnt/lupnt.h>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void init_measurements(py::module &m) {
  m.def(
      "compute_occultation",
      [](double epoch, const Vec3d &r1, const Vec3d &r2, Frame cs1, Frame cs2,
         const std::vector<BodyId> &bodies, const VecXd &atm_h, const VecXd &min_elev,
         const bool use_elev_mask1, const bool use_elev_mask2) -> std::map<std::string, bool> {
        std::map<std::string, bool> vis = Occultation::ComputeOccultation(
            epoch, r1.cast<Real>().eval(), r2.cast<Real>().eval(), cs1, cs2, bodies, atm_h,
            min_elev, use_elev_mask1, use_elev_mask2);
        return vis;
      },
      "Compute occultation between two points", py::arg("epoch"), py::arg("r1"), py::arg("r2"),
      py::arg("cs1"), py::arg("cs2"), py::arg("bodies"), py::arg("atm_h"), py::arg("min_elev"),
      py::arg("use_elev_mask1") = false, py::arg("use_elev_mask2") = false);
  m.def(
      "compute_occultation",
      [](double epoch, Matd<-1, 3> &r1, Matd<-1, 3> &r2, Frame cs1, Frame cs2,
         const std::vector<BodyId> &bodies, const VecXd &atm_h, const VecXd &min_elev,
         const bool use_elev_mask1,
         const bool use_elev_mask2) -> std::vector<std::map<std::string, bool>> {
        return Occultation::ComputeOccultation(epoch, r1.cast<Real>().eval(),
                                               r2.cast<Real>().eval(), cs1, cs2, bodies, atm_h,
                                               min_elev, use_elev_mask1, use_elev_mask2);
      },
      "Compute occultation between two points", py::arg("epoch"), py::arg("r1"), py::arg("r2"),
      py::arg("cs1"), py::arg("cs2"), py::arg("bodies"), py::arg("atm_h"), py::arg("min_elev"),
      py::arg("use_elev_mask1") = false, py::arg("use_elev_mask2") = false);
  m.def(
      "compute_occultation",
      [](const VecXd &epoch, const Matd<-1, 3> &r1, const Matd<-1, 3> &r2, Frame cs1, Frame cs2,
         const std::vector<BodyId> &bodies, const VecXd &atm_h, const VecXd &min_elev,
         const bool use_elev_mask1,
         const bool use_elev_mask2) -> std::vector<std::map<std::string, bool>> {
        return Occultation::ComputeOccultation(epoch.cast<Real>().eval(), r1.cast<Real>().eval(),
                                               r2.cast<Real>().eval(), cs1, cs2, bodies, atm_h,
                                               min_elev, use_elev_mask1, use_elev_mask2);
      },
      "Compute occultation between two points", py::arg("epoch"), py::arg("r1"), py::arg("r2"),
      py::arg("cs1"), py::arg("cs2"), py::arg("bodies"), py::arg("atm_h"), py::arg("min_elev"),
      py::arg("use_elev_mask1") = false, py::arg("use_elev_mask2") = false);
}
