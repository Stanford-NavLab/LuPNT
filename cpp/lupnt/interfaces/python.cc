#ifdef LUPNT_WITH_PYTHON
#  include "lupnt/core/python.h"

namespace lupnt {

  py::module_ py_pnt = []() {
    PythonInterpreter::GetInstance();
    return py::module_::import("pylupnt");
  }();
  py::module_ py_plot = []() {
    PythonInterpreter::GetInstance();
    return py::module_::import("pylupnt.plot");
  }();
  py::module_ py_pio = []() {
    PythonInterpreter::GetInstance();
    return py::module_::import("plotly.io");
  }();
  py::module_ py_go = []() {
    PythonInterpreter::GetInstance();
    return py::module_::import("plotly.graph_objects");
  }();
  py::module_ py_plt = []() {
    PythonInterpreter::GetInstance();
    return py::module_::import("matplotlib.pyplot");
  }();

}  // namespace lupnt
#endif  // LUPNT_WITH_PYTHON
