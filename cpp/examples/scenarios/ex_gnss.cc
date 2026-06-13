// Legacy GNSS scenario placeholder.
//
// The old example used the removed SpaceChannel/SpaceReceiver API. The current
// supported GNSS flow is implemented in the Lunar GNSS ODTS simulation example,
// which precomputes light-time-aware links and then runs filtering from the
// generated data products.
#include <iostream>

int main() {
  std::cout << "The legacy GNSS channel example has been retired.\n"
            << "Use ./build/examples/ex_lunar_gnss_odts --help for the current GNSS "
               "simulation entry point.\n";
  return 0;
}
