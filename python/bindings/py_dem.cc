#include "lupnt/lupnt.h"
#include "py_pybind11.h"
namespace py = pybind11;

void InitDem(py::module& m) {
  m.def("load_tiff", &LoadTiff, py::arg("path"), py::arg("xlims"), py::arg("ylims"),
        py::arg("max_res"));
}
