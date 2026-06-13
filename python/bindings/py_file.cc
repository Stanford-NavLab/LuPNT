
#include <lupnt/lupnt.h>

#include <filesystem>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitFile(py::module& m) {
  m.def("get_data_path", []() { return GetDataPath().string(); });
  m.def("get_output_dir", [](std::string output_dir) { return GetOutputDir(output_dir).string(); });
  m.def("find_file_in_dir", [](std::string base_path, std::string filename) {
    return FindFileInDir(base_path, filename).value_or(std::filesystem::path(""));
  });
  m.def("find_file_in_dir", [](std::string base_path, std::string filename) {
    return FindFileInDir(base_path, filename).value_or(std::filesystem::path(""));
  });
  m.def("get_file_path", [](std::string filename) { return GetFilePath(filename).string(); });
  m.def("get_cspice_kernel_dir", []() { return GetCspiceKernelDir().string(); });
  m.def("get_ascii_kernel_dir", []() { return GetAsciiKernelDir().string(); });
}
