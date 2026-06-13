// Minimal Simulation loader. Pass a YAML config path on the command line:
//   ./build/examples/ex_sim path/to/config.yaml
#include <yaml-cpp/yaml.h>

#include "lupnt/core/logger.h"
#include "lupnt/core/simulation.h"

using namespace lupnt;

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <config.yaml>\n";
    return 2;
  }

  YAML::Node config = YAML::LoadFile(argv[1]);
  Simulation sim(config);
  sim.Run();
  return 0;
}
