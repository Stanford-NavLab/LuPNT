#pragma once

#include <lupnt/core/definitions.h>

// Eigen includes
#include <Eigen/Core>

// pybind11 includes
#include "py_pybind11.h"

template <typename Vec, typename T> void ExportVector(py::module& m, const char* typestr) {
  auto cls = py::class_<Vec>(m, typestr);

  using VecRef = Eigen::Ref<Vec>;
  using VecConstRef = Eigen::Ref<const Vec>;

  cls.def(py::init<>());
  cls.def(py::init<long>());
  cls.def(py::init<const VecRef&>());
  cls.def(py::init<const VecConstRef&>());

  cls.def(py::init([](py::iterable seq) {
            Vec v(py::len(seq));
            size_t i = 0;
            for (auto item : seq) v[i++] = item.cast<T>();
            return v;
          }),
          py::return_value_policy::move);

  cls.def(py::init([](py::buffer buf) {
            py::buffer_info info = buf.request();
            Vec v(info.size);
            auto* src = static_cast<const T*>(info.ptr);
            for (ssize_t i = 0; i < info.size; ++i) v[i] = src[i];
            return v;
          }),
          py::return_value_policy::move);

  cls.def("__len__", [](const Vec& s) { return s.size(); });

  cls.def("__getitem__", [](const Vec& s, size_t i) {
    const size_t size = s.size();
    if (i >= size) throw py::index_error();
    return s[i];
  });

  cls.def("__setitem__", [](Vec& s, size_t i, const T& val) {
    const size_t size = s.size();
    if (i >= size) throw py::index_error();
    s[i] = val;
  });

  cls.def(
      "__iter__", [](const Vec& s) { return py::make_iterator(s.begin(), s.end()); },
      py::keep_alive<0, 1>());

  cls.def("__str__", [](const Vec& s) {
    std::stringstream stream;
    stream << s.transpose().format(lupnt::FMT_COMPACT);
    return stream.str();
  });

  cls.def("__repr__", [=](const Vec& s) {
    std::stringstream stream;
    stream << typestr << "(";
    stream << s.transpose().format(lupnt::FMT_HEAVY);
    stream << ")";
    return stream.str();
  });

  cls.def("__neg__", [](const Vec& s) { return (-s).eval(); });

  cls.def("__add__", [](const Vec& l, const Vec& r) { return (l + r).eval(); });
  cls.def("__add__", [](const T& l, const Vec& r) { return (l + r.array()).matrix().eval(); });
  cls.def("__add__", [](const Vec& l, const T& r) { return (l.array() + r).matrix().eval(); });

  cls.def("__sub__", [](const Vec& l, const Vec& r) { return (l - r).eval(); });
  cls.def("__sub__", [](const T& l, const Vec& r) { return (l - r.array()).matrix().eval(); });
  cls.def("__sub__", [](const Vec& l, const T& r) { return (l.array() - r).matrix().eval(); });

  cls.def("__iadd__", [](Vec& l, const Vec& r) { return (l += r).eval(); });
  cls.def("__iadd__", [](Vec& l, const T& r) { return (l.array() += r).matrix().eval(); });

  cls.def("__isub__", [](Vec& l, const Vec& r) { return (l -= r).eval(); });
  cls.def("__isub__", [](Vec& l, const T& r) { return (l.array() -= r).matrix().eval(); });

  cls.def("__imul__",
          [](Vec& l, const Vec& r) { return (l.array() *= r.array()).matrix().eval(); });
  cls.def("__imul__", [](Vec& l, const T& r) { return (l.array() *= r).matrix().eval(); });

  cls.def("__itruediv__",
          [](Vec& l, const Vec& r) { return (l.array() /= r.array()).matrix().eval(); });
  cls.def("__itruediv__", [](Vec& l, const T& r) { return (l.array() /= r).matrix().eval(); });

  cls.def("__eq__", [](const Vec& l, const Vec& r) { return (l.array() == r.array()).all(); });
  cls.def("__eq__", [](const Vec& l, const T& r) { return (l.array() == r).all(); });

  cls.def("__ne__", [](const Vec& l, const Vec& r) { return (l.array() != r.array()).any(); });
  cls.def("__ne__", [](const Vec& l, const T& r) { return (l.array() != r).any(); });

  cls.def("__mul__",
          [](const Vec& l, const Vec& r) { return (l.array() * r.array()).matrix().eval(); });
  cls.def("__mul__", [](const T& l, const Vec& r) { return (l * r.array()).matrix().eval(); });
  cls.def("__mul__", [](const Vec& l, const T& r) { return (l.array() * r).matrix().eval(); });

  cls.def("__truediv__",
          [](const Vec& l, const Vec& r) { return (l.array() / r.array()).matrix().eval(); });
  cls.def("__truediv__", [](const T& l, const Vec& r) { return (l / r.array()).matrix().eval(); });
  cls.def("__truediv__", [](const Vec& l, const T& r) { return (l.array() / r).matrix().eval(); });

  cls.def("__matmul__", [](const Vec& l, const Vec& r) { return (r.dot(l)); });
  cls.def("__matmul__",
          [](const Vec& l, const T& r) { return (l.transpose() * r).matrix().eval(); });
  cls.def("__matmul__", [](const T& l, const Vec& r) { return (l * r.transpose()).eval(); });

  cls.def("__abs__", [](const Vec& s) { return s.array().abs().matrix().eval(); });

  cls.def(
      "__pow__",
      [](const Vec& base, const Vec& exp) -> Vec {
        return (base.array().pow(exp.array())).matrix();
      },
      py::is_operator(), py::return_value_policy::move);

  cls.def(
      "__pow__",
      [](const Vec& base, const T& exp) -> Vec { return (base.array().pow(exp)).matrix(); },
      py::is_operator(), py::return_value_policy::move);

  cls.def(
      "__rpow__",
      [](const Vec& exp, const T& base) -> Vec {
        return (Eigen::Array<T, Eigen::Dynamic, 1>::Constant(exp.size(), base).pow(exp.array()))
            .matrix();
      },
      py::is_operator(), py::return_value_policy::move);

  cls.def("norm", [](const Vec& v) { return v.norm(); });

  cls.def("dot", [](const Vec& l, const Vec& r) {
    if (l.size() != r.size())
      throw std::runtime_error("Vectors must be of the same size for dot product.");
    return l.dot(r);
  });

  cls.def("sum", [](const Vec& v) { return v.sum(); });

  cls.def("mean", [](const Vec& v) { return v.mean(); });

  cls.def("numpy", [](const Vec& v) {
    py::array_t<double> arr(v.size());
    auto* dst = arr.mutable_data();
    for (Eigen::Index i = 0; i < v.size(); ++i) dst[i] = static_cast<double>(v[i]);
    return arr;
  });
}

template <typename Mat, typename T> void ExportMatrix(py::module& m, const char* typestr) {
  auto cls = py::class_<Mat>(m, typestr);

  using MatRef = Eigen::Ref<Mat>;
  using MatConstRef = Eigen::Ref<const Mat>;

  cls.def(py::init<>());
  cls.def(py::init<long, long>());
  cls.def(py::init<const MatRef&>());
  cls.def(py::init<const MatConstRef&>());

  cls.def("__len__", [](const Mat& s) { return s.size(); });

  cls.def("rows", [](const Mat& s) { return s.rows(); });
  cls.def("cols", [](const Mat& s) { return s.cols(); });

  cls.def("__getitem__", [](const Mat& s, py::tuple pos) {
    const size_t rows = s.rows();
    const size_t cols = s.cols();
    const size_t i = pos[0].cast<size_t>();
    const size_t j = pos[1].cast<size_t>();
    if (i >= rows) throw py::index_error();
    if (j >= cols) throw py::index_error();
    return s(i, j);
  });

  cls.def("__setitem__", [](Mat& s, py::tuple pos, const T& val) {
    const size_t rows = s.rows();
    const size_t cols = s.cols();
    const size_t i = pos[0].cast<size_t>();
    const size_t j = pos[1].cast<size_t>();
    if (i >= rows) throw py::index_error();
    if (j >= cols) throw py::index_error();
    s(i, j) = val;
  });

  cls.def("__str__", [](const Mat& s) {
    std::stringstream stream;
    stream << s.format(lupnt::FMT_COMPACT);
    return stream.str();
  });

  cls.def("__repr__", [=](const Mat& s) {
    std::stringstream stream;
    stream << "pnt." << typestr << "(\n";
    stream << s.format(lupnt::FMT_HEAVY);
    stream << "\n)";
    return stream.str();
  });

  cls.def("numpy", [](const Mat& m) {
    py::array_t<double> arr({m.rows(), m.cols()});
    auto* dst = arr.mutable_data();
    for (Eigen::Index j = 0; j < m.cols(); ++j)
      for (Eigen::Index i = 0; i < m.rows(); ++i) *dst++ = static_cast<double>(m(i, j));
    return arr;
  });
}
