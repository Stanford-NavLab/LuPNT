#include <lupnt/lupnt.h>

using namespace lupnt;
int main() {
  File file("ex_load_save.h5", File::Overwrite);

  MatX A = MatX::Random(10, 5);
  std::cout << A << std::endl << std::endl;
  dump(file, "/path/to/A", A.cast<double>());

  MatX B = load<MatXd>(file, "/path/to/A");
  std::cout << B << std::endl;
}
