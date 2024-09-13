
#include <lupnt/lupnt.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <filesystem>

namespace py = pybind11;
using namespace lupnt;

void init_file(py::module &m) {
  py::class_<std::filesystem::path>(m, "Path")
      .def(py::init<std::string>())
      .def("__str__", [](const std::filesystem::path &p) { return p.string(); })
      .def("__repr__", [](const std::filesystem::path &p) { return p.string(); });
  py::implicitly_convertible<std::string, std::filesystem::path>();

  m.def("get_data_path", &GetDataPath, "Get data path");
  m.def("get_output_path", &GetOutputPath, "Get output path", py::arg("output_dir"));
  m.def("find_file_in_dir", &FindFileInDir, "Find file in directory", py::arg("base_path"),
        py::arg("filename"));
  m.def("get_file_path", &GetFilePath, "Get file path", py::arg("filename"));
  m.def("get_cspice_kernel_dir", &GetCspiceKernelDir, "Get cspice kernel directory");
  m.def("get_ascii_kernel_dir", &GetAsciiKernelDir, "Get ascii kernel directory");
}
