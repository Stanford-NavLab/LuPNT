#include <omp.h>

#include <chrono>
#include <iostream>
#include <thread>

#include "lupnt/environment/plasma/plasma.h"

using namespace pecsim;
using namespace std;

int main(int argc, char** argv) {
  // Lightweight OpenMP smoke test for the shared ProgressBar utility. The
  // optional arguments keep the example fast while still exercising parallel
  // updates:
  //   ex_plasma_openmp [steps] [sleep_ms]
  const int steps = argc > 1 ? std::stoi(argv[1]) : 20;
  const int sleep_ms = argc > 2 ? std::stoi(argv[2]) : 50;

  ProgressBar pg;
  pg.start(steps);

  std::cout << "Max threads: " << omp_get_max_threads() << "\n";
  std::cout << "Num processors: " << omp_get_num_procs() << "\n";

// Be explicit about shared variables when using OpenMP in examples.
#pragma omp parallel for default(none) schedule(static) shared(pg, steps, sleep_ms)
  for (int i = 0; i < steps; i++) {
    pg.update(i);
    std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
  }
}
