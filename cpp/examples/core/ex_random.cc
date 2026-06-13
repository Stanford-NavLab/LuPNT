// Demonstrates deterministic random-number generation by resetting the shared
// LuPNT random engine seed.
#include <lupnt/lupnt.h>
using namespace lupnt;

int main() {
  lupnt::RandomEngine::SetSeed(0);
  std::cout << lupnt::SampleNormal(0, 1) << std::endl;
  std::cout << lupnt::SampleNormal(0, 1) << std::endl;

  lupnt::RandomEngine::SetSeed(0);
  std::cout << lupnt::SampleNormal(0, 1) << std::endl;
  std::cout << lupnt::SampleNormal(0, 1) << std::endl;
}
