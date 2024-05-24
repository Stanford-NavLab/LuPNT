#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

using Eigen::MatrixXd;

// Function to save the entire matrix to a binary file at once
void saveMatrixBuffered(const std::string& filename,
                        const std::vector<MatrixXd>& matrices) {
  std::ofstream out(filename, std::ios::binary);
  if (out.is_open()) {
    for (const auto& matrix : matrices) {
      int rows = matrix.rows();
      int cols = matrix.cols();
      out.write(reinterpret_cast<const char*>(&rows), sizeof(int));
      out.write(reinterpret_cast<const char*>(&cols), sizeof(int));
      out.write(reinterpret_cast<const char*>(matrix.data()),
                rows * cols * sizeof(double));
    }
    out.close();
  } else {
    throw std::runtime_error("Could not open file for writing.");
  }
}

// Function to save a matrix incrementally with the file kept open
void saveMatrixIncremental(std::ofstream& out, const MatrixXd& matrix) {
  if (out.is_open()) {
    int rows = matrix.rows();
    int cols = matrix.cols();
    out.write(reinterpret_cast<const char*>(&rows), sizeof(int));
    out.write(reinterpret_cast<const char*>(&cols), sizeof(int));
    out.write(reinterpret_cast<const char*>(matrix.data()),
              rows * cols * sizeof(double));
  } else {
    throw std::runtime_error("File stream is not open.");
  }
}

int main() {
  int n = 100, m = 100, it = 100000;
  std::ofstream out("incremental_kept_open_matrices.bin", std::ios::binary);
  if (!out.is_open()) {
    std::cerr << "Could not open file for writing." << std::endl;
    return 1;
  }

  auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < it; ++i) {
    MatrixXd mat = MatrixXd::Random(n, m);
    saveMatrixIncremental(out, mat);
  }
  out.close();
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Elapsed time for incremental saving: " << elapsed.count()
            << " s" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  std::vector<MatrixXd> matrices;
  for (int i = 0; i < it; ++i) {
    MatrixXd mat = MatrixXd::Random(n, m);
    matrices.push_back(mat);
  }
  saveMatrixBuffered("buffered_matrices.bin", matrices);
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  std::cout << "Elapsed time for buffered saving: " << elapsed.count() << " s"
            << std::endl;

  return 0;
}
