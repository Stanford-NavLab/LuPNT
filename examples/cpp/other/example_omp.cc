#include <omp.h>

#include <iostream>

int main() {
#pragma omp parallel for
  for (int i = 0; i < 10; i++) {
    int thread_id = omp_get_thread_num();
    int num_threads = omp_get_num_threads();
    std::cout << "Hello from thread " << thread_id << " out of " << num_threads << " threads.\n";
  }
  return 0;
}
