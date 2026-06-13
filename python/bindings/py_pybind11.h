#pragma once

//----------------------------------------------------------------------------------------------------------------
// Add defines above to allow pybind11 packages produced with different compilers/versions to be
// used together. https://github.com/pybind/pybind11/pull/2602
//----------------------------------------------------------------------------------------------------------------
// #define PYBIND11_COMPILER_TYPE ""
// #define PYBIND11_STDLIB ""
// #define PYBIND11_BUILD_ABI ""
//----------------------------------------------------------------------------------------------------------------

// Include lupnt definitions first to get Real type
#include <lupnt/core/definitions.h>
#include <lupnt/core/error.h>
using namespace lupnt;

// Include basic pybind11 headers first
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

// Include OpenCV type caster
// #include "py_opencv.h"  // NOLINT(misc-include-cleaner) - provides cv::Mat type caster

namespace pybind11 {
  namespace detail {
    template <> struct type_caster<Real> {
      PYBIND11_TYPE_CASTER(Real, const_name("Real"));
      bool load(handle src, bool) {
        // Only handle scalar conversions, not arrays
        if (isinstance<array>(src)) {
          return false;  // Let Eigen handle array conversions
        }
        PyObject* tmp = PyNumber_Float(src.ptr());
        if (!tmp) return false;
        double v = PyFloat_AsDouble(tmp);
        Py_DECREF(tmp);
        if (v == -1.0 && PyErr_Occurred()) return false;
        value = Real(v);
        return true;
      }
      static handle cast(const Real& src, return_value_policy, handle) {
        return PyFloat_FromDouble(static_cast<double>(src));
      }
    };

  }  // namespace detail
}  // namespace pybind11

// Define npy_format_descriptor for Real before including eigen
namespace pybind11 {
  namespace detail {
    template <> struct npy_format_descriptor<Real> {
      static constexpr auto name = _("Real");
      static std::string format() { return format_descriptor<double>::format(); }
      static pybind11::dtype dtype() { return pybind11::dtype::of<double>(); }
    };
  }  // namespace detail
}  // namespace pybind11

// Now define Eigen dtype and include eigen
#define PYBIND11_EIGEN_DTYPE Real
#include <pybind11/eigen.h>

// Include remaining pybind11 headers
#include <pybind11/operators.h>
#include <pybind11/stl.h>

// Set up namespace aliases
namespace py = pybind11;
