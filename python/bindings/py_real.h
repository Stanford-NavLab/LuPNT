#pragma once

// pybind11 includes
#include "py_pybind11.h"

// C++ includes
#include <sstream>

template <size_t N, typename T> void ExportReal(py::module& m, const char* typestr) {
  auto __getitem__ = [](const ad::Real<N, T>& self, size_t i) { return self[i]; };

  auto __setitem__ = [](ad::Real<N, T>& self, size_t i, const T& value) { self[i] = value; };

  auto __str__ = [](const ad::Real<N, T>& self) {
    std::stringstream ss;
    ss << self;
    return ss.str();
  };

  auto __repr__ = [](const ad::Real<N, T>& self) {
    std::stringstream ss;
    ss << "Real(";
    for (size_t i = 0; i < N; ++i) ss << self[i] << ", ";
    ss << self[N] << ")";
    return ss.str();
  };

  auto __float__ = [](const ad::Real<N, T>& self) { return self[0]; };

  auto __deepcopy__
      = [](const ad::Real<N, T>& self, py::object memo) -> ad::Real<N, T> { return self; };

  auto cls
      = py::class_<ad::Real<N, T>>(m, typestr)
            .def(py::init<>())
            .def(py::init<const T&>())
            .def(py::init<const std::array<T, N + 1>&>())
            .def(py::init<const ad::Real<N, T>&>())
            .def("__getitem__", __getitem__)
            .def("__setitem__", __setitem__)
            .def("__str__", __str__)
            .def("__repr__", __repr__)
            .def("__float__", __float__)
            .def("__deepcopy__",
                 __deepcopy__)  // needed when using autodiff::real types with plotly

            .def("exp", &autodiff::detail::exp<N, T>)
            .def("log", &autodiff::detail::log<N, T>)
            .def("log10", &autodiff::detail::log10<N, T>)
            .def("sqrt", &autodiff::detail::sqrt<N, T>)
            .def("cbrt", &autodiff::detail::cbrt<N, T>)
            .def("__pow__", &autodiff::detail::pow<N, T>)
            .def("__pow__", [](const ad::Real<N, T>& self, const T& x) { return pow(self, x); })
            .def("__rpow__", [](const ad::Real<N, T>& self, const T& x) { return pow(x, self); })
            .def("sin", &autodiff::detail::sin<N, T>)
            .def("cos", &autodiff::detail::cos<N, T>)
            .def("tan", &autodiff::detail::tan<N, T>)
            .def("arcsin", &autodiff::detail::asin<N, T>)
            .def("arccos", &autodiff::detail::acos<N, T>)
            .def("arctan", &autodiff::detail::atan<N, T>)
            .def("arctan2", &autodiff::detail::atan2<N, T>)
            .def("arctan2", [](const ad::Real<N, T>& x, const T& y) { return atan2(x, y); })
            .def("sinh", &autodiff::detail::sinh<N, T>)
            .def("cosh", &autodiff::detail::cosh<N, T>)
            .def("tanh", &autodiff::detail::tanh<N, T>)
            .def("arcsinh", &autodiff::detail::asinh<N, T>)
            .def("arccosh", &autodiff::detail::acosh<N, T>)
            .def("__abs__", &autodiff::detail::abs<N, T>)
            .def("minimum", &autodiff::detail::min<N, T>)
            .def("minimum", [](const ad::Real<N, T>& x, const T& y) { return min(x, y); })
            .def("minimum", [](const T& x, const ad::Real<N, T>& y) { return min(x, y); })
            .def("maximum", &autodiff::detail::max<N, T>)
            .def("maximum", [](const ad::Real<N, T>& x, const T& y) { return max(x, y); })
            .def("maximum", [](const T& x, const ad::Real<N, T>& y) { return max(x, y); })

            .def("val", [](const ad::Real<N, T>& self) { return self.val(); })
            .def("val", [](ad::Real<N, T>& self) { return self.val(); })

            .def(-py::self)

            .def(py::self + py::self)
            .def(py::self - py::self)
            .def(py::self * py::self)
            .def(py::self / py::self)

            .def(py::self + T())
            .def(py::self - T())
            .def(py::self * T())
            .def(py::self / T())

            .def(T() + py::self)
            .def(T() - py::self)
            .def(T() * py::self)
            .def(T() / py::self)

            .def(py::self += py::self)
            .def(py::self -= py::self)
            .def(py::self *= py::self)
            .def(py::self /= py::self)

            .def(py::self == py::self)
            .def(py::self != py::self)
            .def(py::self < py::self)
            .def(py::self > py::self)
            .def(py::self <= py::self)
            .def(py::self >= py::self)

            .def(py::self += T())
            .def(py::self -= T())
            .def(py::self *= T())
            .def(py::self /= T())

            .def(py::self == T())
            .def(py::self != T())
            .def(py::self < T())
            .def(py::self > T())
            .def(py::self <= T())
            .def(py::self >= T())

            .def(T() == py::self)
            .def(T() != py::self)
            .def(T() < py::self)
            .def(T() > py::self)
            .def(T() <= py::self)
            .def(T() >= py::self);

  if constexpr (!isSame<T, int>) cls.def(py::init<int>());
  if constexpr (!isSame<T, long>) cls.def(py::init<long>());
  if constexpr (!isSame<T, float>) cls.def(py::init<float>());
  if constexpr (!isSame<T, double>) cls.def(py::init<double>());

  py::implicitly_convertible<T, ad::Real<N, T>>();

  if constexpr (!isSame<T, int>) py::implicitly_convertible<int, ad::Real<N, T>>();
  if constexpr (!isSame<T, long>) py::implicitly_convertible<long, ad::Real<N, T>>();
  if constexpr (!isSame<T, float>) py::implicitly_convertible<float, ad::Real<N, T>>();
  if constexpr (!isSame<T, double>) py::implicitly_convertible<double, ad::Real<N, T>>();

  m.def("seed", [](ad::Real<N, T>& x) { x[1] = 1.0; });
  m.def("unseed", [](ad::Real<N, T>& x) { x[1] = 1.0; });

  m.def("abs", [](const ad::Real<N, T>& x) { return abs(x); });

  m.def("sin", [](const ad::Real<N, T>& x) { return sin(x); });
  m.def("cos", [](const ad::Real<N, T>& x) { return cos(x); });
  m.def("tan", [](const ad::Real<N, T>& x) { return tan(x); });

  m.def("asin", [](const ad::Real<N, T>& x) { return asin(x); });
  m.def("acos", [](const ad::Real<N, T>& x) { return acos(x); });
  m.def("atan", [](const ad::Real<N, T>& x) { return atan(x); });

  m.def("asinh", [](const ad::Real<N, T>& x) { return asinh(x); });
  m.def("acosh", [](const ad::Real<N, T>& x) { return acosh(x); });
  m.def("atanh", [](const ad::Real<N, T>& x) { return atanh(x); });

  m.def("sinh", [](const ad::Real<N, T>& x) { return sinh(x); });
  m.def("cosh", [](const ad::Real<N, T>& x) { return cosh(x); });
  m.def("tanh", [](const ad::Real<N, T>& x) { return tanh(x); });

  m.def("arcsin", [](const ad::Real<N, T>& x) { return asin(x); });
  m.def("arccos", [](const ad::Real<N, T>& x) { return acos(x); });
  m.def("arctan", [](const ad::Real<N, T>& x) { return atan(x); });

  m.def("arcsinh", [](const ad::Real<N, T>& x) { return asinh(x); });
  m.def("arccosh", [](const ad::Real<N, T>& x) { return acosh(x); });
  m.def("arctanh", [](const ad::Real<N, T>& x) { return atanh(x); });

  m.def("sqrt", [](const ad::Real<N, T>& x) { return sqrt(x); });
  m.def("cbrt", [](const ad::Real<N, T>& x) { return cbrt(x); });

  m.def("exp", [](const ad::Real<N, T>& x) { return exp(x); });
  m.def("log", [](const ad::Real<N, T>& x) { return log(x); });
  m.def("log10", [](const ad::Real<N, T>& x) { return log10(x); });

  m.def("pow", [](const ad::Real<N, T>& x, const ad::Real<N, T>& y) { return pow(x, y); });
  m.def("pow", [](const ad::Real<N, T>& x, const T& y) { return pow(x, y); });
  m.def("pow", [](const T& x, const ad::Real<N, T>& y) { return pow(x, y); });

  m.def("max", [](const ad::Real<N, T>& x, const ad::Real<N, T>& y) { return max(x, y); });
  m.def("max", [](const ad::Real<N, T>& x, const T& y) { return max(x, y); });
  m.def("max", [](const T& x, const ad::Real<N, T>& y) { return max(x, y); });

  m.def("min", [](const ad::Real<N, T>& x, const ad::Real<N, T>& y) { return min(x, y); });
  m.def("min", [](const ad::Real<N, T>& x, const T& y) { return min(x, y); });
  m.def("min", [](const T& x, const ad::Real<N, T>& y) { return min(x, y); });
}
