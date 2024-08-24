#define VEC_BIND_REAL(name, func, arg1)                                                         \
  m.def(name, [](double x) -> double { return func(x).val(); }, py::arg(arg1));                 \
  m.def(                                                                                        \
      name, [](const VecXd &x) -> VecXd { return func(x.cast<Real>().eval()).cast<double>(); }, \
      py::arg(arg1));

#define VEC_BIND_VECTOR(name, func, size, arg1)            \
  m.def(                                                   \
      name,                                                \
      [](const Vecd<size> &x) -> Vecd<size> {              \
        return func(x.cast<Real>().eval()).cast<double>(); \
      },                                                   \
      py::arg(arg1));                                      \
  m.def(                                                   \
      name,                                                \
      [](const Matd<-1, size> &x) -> Matd<-1, size> {      \
        return func(x.cast<Real>().eval()).cast<double>(); \
      },                                                   \
      py::arg(arg1));

#define VEC_BIND_VECTOR_REAL(name, func, size, arg1, arg2)                        \
  m.def(                                                                          \
      name,                                                                       \
      [](const Vecd<size> &x, double y) -> Vecd<size> {                           \
        return func(x.cast<Real>().eval(), y).cast<double>();                     \
      },                                                                          \
      py::arg(arg1), py::arg(arg2));                                              \
  m.def(                                                                          \
      name,                                                                       \
      [](const Vecd<size> &x, const VecXd &y) -> Matd<-1, size> {                 \
        return func(x.cast<Real>().eval(), y.cast<Real>().eval()).cast<double>(); \
      },                                                                          \
      py::arg(arg1), py::arg(arg2));                                              \
  m.def(                                                                          \
      name,                                                                       \
      [](const Matd<-1, size> &x, double y) -> Matd<-1, size> {                   \
        return func(x.cast<Real>().eval(), y).cast<double>();                     \
      },                                                                          \
      py::arg(arg1), py::arg(arg2));                                              \
  m.def(                                                                          \
      name,                                                                       \
      [](const Matd<-1, size> &x, const VecXd &y) -> Matd<-1, size> {             \
        return func(x.cast<Real>().eval(), y.cast<Real>().eval()).cast<double>(); \
      },                                                                          \
      py::arg(arg1), py::arg(arg2));

#define VEC_BIND_VECTOR_VECTOR(name, func, size, arg1, arg2)                              \
  m.def(                                                                                  \
      name,                                                                               \
      [](const Vecd<size> &x, const Vecd<size> &y) -> Vecd<size> {                        \
        return func(x.cast<Real>().eval(), y.cast<Real>().eval()).cast<double>();         \
      },                                                                                  \
      py::arg(arg1), py::arg(arg2));                                                      \
  m.def(                                                                                  \
      name,                                                                               \
      [](const RowVecd<size> &x, const RowVecd<size> &y) -> RowVecd<size> {               \
        return func(x.transpose().cast<Real>().eval(), y.transpose().cast<Real>().eval()) \
            .cast<double>();                                                              \
      },                                                                                  \
      py::arg(arg1), py::arg(arg2));                                                      \
  m.def(                                                                                  \
      name,                                                                               \
      [](const Matd<-1, size> &x, const Vecd<size> &y) -> Matd<-1, size> {                \
        return func(x.cast<Real>().eval(), y.cast<Real>().eval()).cast<double>();         \
      },                                                                                  \
      py::arg(arg1), py::arg(arg2));                                                      \
  m.def(                                                                                  \
      name,                                                                               \
      [](const Vecd<size> &x, const Matd<-1, size> &y) -> Matd<-1, size> {                \
        return func(x.cast<Real>().eval(), y.cast<Real>().eval()).cast<double>();         \
      },                                                                                  \
      py::arg(arg1), py::arg(arg2));                                                      \
  m.def(                                                                                  \
      name,                                                                               \
      [](const Matd<-1, size> &x, const Matd<-1, size> &y) -> Matd<-1, size> {            \
        return func(x.cast<Real>().eval(), y.cast<Real>().eval()).cast<double>();         \
      },                                                                                  \
      py::arg(arg1), py::arg(arg2));

#define VEC_BIND_REAL_REAL(name, func, arg1, arg2)                                             \
  m.def(                                                                                       \
      name, [](double x, double y) -> double { return func(x, y).val(); }, py::arg(arg1),      \
      py::arg(arg2));                                                                          \
  m.def(                                                                                       \
      name,                                                                                    \
      [](const VecXd &x, double y) -> VecXd {                                                  \
        return func(x.cast<Real>().eval(), y).cast<double>();                                  \
      },                                                                                       \
      py::arg(arg1), py::arg(arg2));                                                           \
  m.def(                                                                                       \
      name, [](double x, const VecXd &y) -> VecXd { return func(x, y).cast<double>(); },       \
      py::arg(arg1), py::arg(arg2));                                                           \
  m.def(                                                                                       \
      name, [](const VecXd &x, const VecXd &y) -> VecXd { return func(x, y).cast<double>(); }, \
      py::arg(arg1), py::arg(arg2));
