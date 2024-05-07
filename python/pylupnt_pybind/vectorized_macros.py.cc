#define VECTORIZED_BINDING_FROM_VECTOR(name, func, size, arg1) \
  m.def(                                                       \
      name,                                                    \
      [](const Vectord<size> &x) -> Vectord<size> {            \
        return func(x.cast<real>().eval()).cast<double>();     \
      },                                                       \
      py::arg(arg1));                                          \
  m.def(                                                       \
      name,                                                    \
      [](const Matrixd<-1, size> &x) -> Matrixd<-1, size> {    \
        return func(x.cast<real>().eval()).cast<double>();     \
      },                                                       \
      py::arg(arg1));

#define VECTORIZED_BINDING_FROM_VECTOR_REAL(name, func, size, arg1, arg2)      \
  m.def(                                                                       \
      name,                                                                    \
      [](const Vectord<size> &x, double y) -> Vectord<size> {                  \
        return func(x.cast<real>().eval(), y).cast<double>();                  \
      },                                                                       \
      py::arg(arg1), py::arg(arg2));                                           \
  m.def(                                                                       \
      name,                                                                    \
      [](const Vectord<size> &x, const VectorXd &y) -> Matrixd<-1, size> {     \
        return func(x.cast<real>().eval(), y.cast<real>().eval())              \
            .cast<double>();                                                   \
      },                                                                       \
      py::arg(arg1), py::arg(arg2));                                           \
  m.def(                                                                       \
      name,                                                                    \
      [](const Matrixd<-1, size> &x, double y) -> Matrixd<-1, size> {          \
        return func(x.cast<real>().eval(), y).cast<double>();                  \
      },                                                                       \
      py::arg(arg1), py::arg(arg2));                                           \
  m.def(                                                                       \
      name,                                                                    \
      [](const Matrixd<-1, size> &x, const VectorXd &y) -> Matrixd<-1, size> { \
        return func(x.cast<real>().eval(), y.cast<real>().eval())              \
            .cast<double>();                                                   \
      },                                                                       \
      py::arg(arg1), py::arg(arg2));

#define VECTORIZED_BINDING_FROM_VECTOR_VECTOR(name, func, size, arg1, arg2) \
  m.def(                                                                    \
      name,                                                                 \
      [](const Vectord<size> &x, const Vectord<size> &y) -> Vectord<size> { \
        return func(x.cast<real>().eval(), y.cast<real>().eval())           \
            .cast<double>();                                                \
      },                                                                    \
      py::arg(arg1), py::arg(arg2));                                        \
  m.def(                                                                    \
      name,                                                                 \
      [](const RowVectord<size> &x,                                         \
         const RowVectord<size> &y) -> RowVectord<size> {                   \
        return func(x.transpose().cast<real>().eval(),                      \
                    y.transpose().cast<real>().eval())                      \
            .cast<double>();                                                \
      },                                                                    \
      py::arg(arg1), py::arg(arg2));                                        \
  m.def(                                                                    \
      name,                                                                 \
      [](const Matrixd<-1, size> &x,                                        \
         const Vectord<size> &y) -> Matrixd<-1, size> {                     \
        return func(x.cast<real>().eval(), y.cast<real>().eval())           \
            .cast<double>();                                                \
      },                                                                    \
      py::arg(arg1), py::arg(arg2));                                        \
  m.def(                                                                    \
      name,                                                                 \
      [](const Vectord<size> &x,                                            \
         const Matrixd<-1, size> &y) -> Matrixd<-1, size> {                 \
        return func(x.cast<real>().eval(), y.cast<real>().eval())           \
            .cast<double>();                                                \
      },                                                                    \
      py::arg(arg1), py::arg(arg2));                                        \
  m.def(                                                                    \
      name,                                                                 \
      [](const Matrixd<-1, size> &x,                                        \
         const Matrixd<-1, size> &y) -> Matrixd<-1, size> {                 \
        return func(x.cast<real>().eval(), y.cast<real>().eval())           \
            .cast<double>();                                                \
      },                                                                    \
      py::arg(arg1), py::arg(arg2));

#define VECTORIZED_BINDING_FROM_REAL_REAL(name, func, arg1, arg2)          \
  m.def(                                                                   \
      name, [](double x, double y) -> double { return func(x, y).val(); }, \
      py::arg(arg1), py::arg(arg2));                                       \
  m.def(                                                                   \
      name,                                                                \
      [](const VectorXd &x, double y) -> VectorXd {                        \
        return func(x.cast<real>().eval(), y).cast<double>();              \
      },                                                                   \
      py::arg(arg1), py::arg(arg2));                                       \
  m.def(                                                                   \
      name,                                                                \
      [](double x, const VectorXd &y) -> VectorXd {                        \
        return func(x, y).cast<double>();                                  \
      },                                                                   \
      py::arg(arg1), py::arg(arg2));                                       \
  m.def(                                                                   \
      name,                                                                \
      [](const VectorXd &x, const VectorXd &y) -> VectorXd {               \
        return func(x, y).cast<double>();                                  \
      },                                                                   \
      py::arg(arg1), py::arg(arg2));