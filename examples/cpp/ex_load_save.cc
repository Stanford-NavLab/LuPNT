#include <lupnt/lupnt.h>

#include <highfive/H5Easy.hpp>

using namespace lupnt;
int main() {
  H5Easy::File file("ex_load_save.h5", H5Easy::File::Overwrite);

  MatX A = MatX::Random(10, 5);
  std::cout << A << std::endl << std::endl;
  H5Easy::dump(file, "/path/to/A", A.cast<double>());

  MatX B = H5Easy::load<MatXd>(file, "/path/to/A");
  std::cout << B << std::endl;
}
