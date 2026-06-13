// OpenMP microbenchmark for a simple vector multiply. Pass an optional element
// count when you want a heavier run:
//   ./build/examples/ex_omp 50000000
#include <omp.h>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
  const std::size_t N_elements = argc > 1 ? std::stoull(argv[1]) : 5'000'000;
  std::vector<double> data(N_elements, 1.0);
  std::vector<double> result(N_elements, 0.0);

  // Serial execution
  auto start_serial = std::chrono::high_resolution_clock::now();
  for (std::size_t i = 0; i < N_elements; ++i) {
    result[i] = data[i] * 2.0;
  }
  auto end_serial = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration_serial = end_serial - start_serial;
  std::cout << "Serial execution time: " << duration_serial.count() << " seconds" << std::endl;

  // Parallel execution with OpenMP
  auto start_parallel = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
  for (std::int64_t i = 0; i < static_cast<std::int64_t>(N_elements); ++i) {
    if (i == 0) {
      int num_threads = omp_get_num_threads();
      std::cout << "Number of threads: " << num_threads << std::endl;
    }
    result[i] = data[i] * 2.0;
  }
  auto end_parallel = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration_parallel = end_parallel - start_parallel;
  std::cout << "Parallel execution time: " << duration_parallel.count() << " seconds" << std::endl;

  // Prevent compiler optimizing away the loop
  double sum = 0.0;
  for (std::size_t i = 0; i < N_elements; ++i) {
    sum += result[i];
  }
  std::cout << "Sum of result: " << sum << std::endl;

  return 0;
}
