
#include <lupnt/lupnt.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <filesystem>

namespace py = pybind11;
using namespace lupnt;

void init_file(py::module& m) {
  m.def("get_data_path", []() { return GetDataPath().string(); });
  m.def("get_output_path",
        [](std::string output_dir) { return GetOutputPath(output_dir).string(); });
  m.def("find_file_in_dir", [](std::string base_path, std::string filename) {
    return FindFileInDir(base_path, filename).value_or(std::filesystem::path(""));
  });
  m.def("get_file_path", [](std::string filename) { return GetFilePath(filename).string(); });
  m.def("get_cspice_kernel_dir", []() { return GetCspiceKernelDir().string(); });
  m.def("get_ascii_kernel_dir", []() { return GetAsciiKernelDir().string(); });
}
