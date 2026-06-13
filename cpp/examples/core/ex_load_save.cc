// Demonstrates HDF5 matrix I/O with HighFive/H5Easy, including attaching and
// reading simple string attributes.
#include <lupnt/lupnt.h>

#include <highfive/H5Easy.hpp>

using namespace lupnt;
int main() {
  H5Easy::File file("ex_load_save.h5", H5Easy::File::Overwrite);

  MatX A = MatX::Random(10, 5);
  std::cout << A << std::endl << std::endl;
  H5Easy::dump(file, "path/to/A", A.cast<double>().eval());
  std::vector<std::string> units = {"km", "m", "s", "rad", "rad/s", "rad/s^2"};
  H5Easy::dumpAttribute(file, "path/to/A", "units", units);

  MatX B = H5Easy::load<MatXd>(file, "path/to/A");
  std::cout << B << std::endl;

  std::cout << H5Easy::loadAttribute<std::vector<std::string>>(file, "path/to/A", "units")
            << std::endl;
}
