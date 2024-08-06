#include <lupnt/lupnt.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

void init_datasets(py::module& m) {
  py::class_<CraterDataLoader>(m, "CraterDataLoader")
      .def_static(
          "load_craters",
          [](const std::string& filename, std::optional<std::pair<double, double>> latlims,
             std::optional<std::pair<double, double>> longlims, std::pair<double, double> diamlims,
             double ellipse_limit, double arc_lims) {
            const std::pair<double, double>* latlims_ptr = latlims ? &*latlims : nullptr;
            const std::pair<double, double>* longlims_ptr = longlims ? &*longlims : nullptr;
            return CraterDataLoader::LoadCraters(filename, latlims_ptr, longlims_ptr, diamlims,
                                                 ellipse_limit, arc_lims);
          },
          py::arg("filename") = "lunar_crater_database_robbins_2018.csv",
          py::arg("latlims") = py::none(), py::arg("longlims") = py::none(),
          py::arg("diamlims") = std::make_pair(0.0, 500.0), py::arg("ellipse_limit") = 1.5,
          py::arg("arc_lims") = 0.0, "Load craters from the dataset")
      .def_static(
          "extract_robbins_dataset",
          [](const std::vector<Crater>& craters, bool radians) {
            Eigen::VectorXd lat, lon, major, minor, psi;
            std::vector<std::string> crater_id;
            CraterDataLoader::ExtractRobbinsDataset(craters, lat, lon, major, minor, psi, crater_id,
                                                    radians);
            return py::make_tuple(lat, lon, major, minor, psi, crater_id);
          },
          py::arg("craters"), py::arg("radians") = true, "Extract the Robbins dataset");

  py::class_<Crater>(m, "Crater")
      .def_readwrite("lat", &Crater::lat)
      .def_readwrite("lon", &Crater::lon)
      .def_readwrite("diam_major", &Crater::diam_major)
      .def_readwrite("diam_minor", &Crater::diam_minor)
      .def_readwrite("angle", &Crater::angle)
      .def_readwrite("id", &Crater::id);
}