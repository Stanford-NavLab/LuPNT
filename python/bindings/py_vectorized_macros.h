#pragma once

#define DEF_REAL(name, func, arg1) def(name, py::overload_cast<Real>(&func), py::arg(arg1))

#define DEF_VECTOR(name, func, size, arg1)                              \
  def(name, py::overload_cast<const Vec<size> &>(&func), py::arg(arg1)) \
      .def(name, py::overload_cast<const Mat<-1, size> &>(&func), py::arg(arg1))

#define DEF_VECTOR_REAL(name, func, size, arg1, arg2)                                          \
  def(name, py::overload_cast<const Vec<size> &, Real>(&func), py::arg(arg1), py::arg(arg2))   \
      .def(name, py::overload_cast<const Vec<size> &, const VecX &>(&func), py::arg(arg1),     \
           py::arg(arg2))                                                                      \
      .def(name, py::overload_cast<const Mat<-1, size> &, Real>(&func), py::arg(arg1),         \
           py::arg(arg2))                                                                      \
      .def(name, py::overload_cast<const Mat<-1, size> &, const VecX &>(&func), py::arg(arg1), \
           py::arg(arg2))

#define DEF_VECTOR_VECTOR(name, func, size, arg1, arg2)                                     \
  def(                                                                                      \
      name, [](const Vec<size> &x1, Vec<size> &x2) { return func(x1, x2); }, py::arg(arg1), \
      py::arg(arg2))                                                                        \
      .def(                                                                                 \
          name, [](const Mat<-1, size> &x1, Vec<size> &x2) { return func(x1, x2); },        \
          py::arg(arg1), py::arg(arg2))                                                     \
      .def(                                                                                 \
          name, [](const Vec<size> &x1, Mat<-1, size> &x2) { return func(x1, x2); },        \
          py::arg(arg1), py::arg(arg2))                                                     \
      .def(                                                                                 \
          name, [](const Mat<-1, size> &x1, Mat<-1, size> &x2) { return func(x1, x2); },    \
          py::arg(arg1), py::arg(arg2))

#define DEF_REAL_REAL(name, func, arg1, arg2)                                                \
  def(name, py::overload_cast<Real, Real>(&func), py::arg(arg1), py::arg(arg2))              \
      .def(name, py::overload_cast<const VecX &, Real>(&func), py::arg(arg1), py::arg(arg2)) \
      .def(name, py::overload_cast<Real, const VecX &>(&func), py::arg(arg1), py::arg(arg2)) \
      .def(name, py::overload_cast<const VecX &, const VecX &>(&func), py::arg(arg1),        \
           py::arg(arg2))

#define DEF_CLASS_REAL_REAL(name, class, func, arg1, arg2)                                         \
  def(                                                                                             \
      name, [](class &cl, Real x, Real y) { return cl.func(x, y); }, py::arg(arg1), py::arg(arg2)) \
      .def(                                                                                        \
          name, [](class &cl, const VecX &x, Real y) { return cl.func(x, y); }, py::arg(arg1),     \
          py::arg(arg2))                                                                           \
      .def(                                                                                        \
          name, [](class &cl, Real x, const VecX &y) { return cl.func(x, y); }, py::arg(arg1),     \
          py::arg(arg2))                                                                           \
      .def(                                                                                        \
          name, [](class &cl, const VecX &x, const VecX &y) { return cl.func(x, y); },             \
          py::arg(arg1), py::arg(arg2))
