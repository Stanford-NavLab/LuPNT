
#include <lupnt/lupnt.h>

#include <filesystem>

#include "py_pybind11.h"

namespace py = pybind11;
using namespace lupnt;

void InitFile(py::module& m) {
  m.def(
      "get_data_path", []() { return GetDataPath().string(); },
      "Path to the LuPNT data directory (LUPNT_DATA_PATH).");
  m.def(
      "get_output_dir", [](std::string output_dir) { return GetOutputDir(output_dir).string(); },
      "Path to a named output subdirectory, creating it if needed.");
  m.def(
      "find_file_in_dir",
      [](std::string base_path, std::string filename) {
        return FindFileInDir(base_path, filename).value_or(std::filesystem::path(""));
      },
      "Recursively search base_path for filename; empty path if not found.");
  m.def(
      "find_file_in_dir",
      [](std::string base_path, std::string filename) {
        return FindFileInDir(base_path, filename).value_or(std::filesystem::path(""));
      },
      "Recursively search base_path for filename; empty path if not found.");
  m.def(
      "get_file_path", [](std::string filename) { return GetFilePath(filename).string(); },
      "Resolve a filename to a full path within the LuPNT data search dirs.");
  m.def(
      "get_cspice_kernel_dir", []() { return GetCspiceKernelDir().string(); },
      "Path to the CSPICE binary kernel directory.");
  m.def(
      "get_ascii_kernel_dir", []() { return GetAsciiKernelDir().string(); },
      "Path to the ASCII (text) SPICE kernel directory.");
}
