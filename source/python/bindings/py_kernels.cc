// lupnt
#include <lupnt/core/constants.h>
#include <lupnt/data/kernels.h>
#include <lupnt/numerics/math_utils.h>

// pybind11
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>

#include "py_vectorized_macros.cc"

namespace py = pybind11;
using namespace lupnt;

void init_kernels(py::module &m) {
  // GetBodyPosVel
  m.def(
      "get_body_pos_vel",
      [](double t_tai, NaifId center, NaifId target, Frame frame) -> Vec6d {
        Real t_tai_ = t_tai;
        return GetBodyPosVel(t_tai_, center, target, frame).cast<double>();
      },
      py::arg("t_tai"), py::arg("center"), py::arg("target"), py::arg("frame"));

  m.def(
      "get_body_pos_vel",
      [](VecXd t_tai, NaifId center, NaifId target, Frame frame) -> MatX6d {
        VecX t_tai_ = t_tai.cast<Real>().array();
        return GetBodyPosVel(t_tai_, center, target, frame).cast<double>();
      },
      py::arg("t_tai"), py::arg("center"), py::arg("target"), py::arg("frame"));

  m.def(
      "get_body_pos_vel",
      [](double t_tai, NaifId target, Frame frame) -> Vec6d {
        Real t_tai_ = t_tai;
        return GetBodyPosVel(t_tai_, target, frame).cast<double>();
      },
      py::arg("t_tai"), py::arg("target"), py::arg("frame"));

  m.def(
      "get_body_pos_vel",
      [](VecXd t_tai, NaifId target, Frame frame) -> MatX6d {
        VecX t_tai_ = t_tai.cast<Real>().array();
        return GetBodyPosVel(t_tai_, target, frame).cast<double>();
      },
      py::arg("t_tai"), py::arg("target"), py::arg("frame"));
}
