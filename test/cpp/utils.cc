
#include <lupnt/core/constants.h>
#include <lupnt/numerics/math_utils.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

using namespace Catch::Matchers;
using namespace lupnt;

#define REQUIRE_NEAR_REAL_VEC(a, b, abs_error)                  \
  REQUIRE(a.size() == b.size());                                \
  for (int i = 0; i < a.size(); ++i) {                          \
    REQUIRE_THAT(a[i].val(), WithinAbs(b[i].val(), abs_error)); \
  }

#define REQUIRE_NEAR_REAL_MAT(a, b, abs_error)                          \
  REQUIRE(a.rows() == b.rows());                                        \
  REQUIRE(a.cols() == b.cols());                                        \
  for (int i = 0; i < a.rows(); ++i) {                                  \
    for (int j = 0; j < a.cols(); ++j) {                                \
      REQUIRE_THAT(a(i, j).val(), WithinAbs(b(i, j).val(), abs_error)); \
    }                                                                   \
  }

#define REQUIRE_NEAR_DOUBLE_VEC(a, b, abs_error)    \
  REQUIRE(a.size() == b.size());                    \
  for (int i = 0; i < a.size(); ++i) {              \
    REQUIRE_THAT(a[i], WithinAbs(b[i], abs_error)); \
  }

#define REQUIRE_NEAR_DOUBLE_MAT(a, b, abs_error)            \
  REQUIRE(a.rows() == b.rows());                            \
  REQUIRE(a.cols() == b.cols());                            \
  for (int i = 0; i < a.rows(); ++i) {                      \
    for (int j = 0; j < a.cols(); ++j) {                    \
      REQUIRE_THAT(a(i, j), WithinAbs(b(i, j), abs_error)); \
    }                                                       \
  }

#define REQUIRE_NEAR_REAL(a, b, abs_error) REQUIRE_THAT(a.val(), WithinAbs(b.val(), abs_error))

#define REQUIRE_NEAR_DOUBLE(a, b, abs_error) REQUIRE_THAT(a, WithinAbs(b, abs_error))

// inputs: vector, function(vector), jacobian
static void NumericalJacobian(std::function<void(VecX&, Real)> propagate_function, const VecX& vec,
                              Real dt, Mat6d& jacobian, double eps = 1e-6) {
  int n = vec.size();
  VecX vec_p;
  VecX vec_m;

  for (int i = 0; i < n; i++) {
    vec_p = vec;
    vec_m = vec;
    vec_p(i) += eps;
    vec_m(i) -= eps;

    propagate_function(vec_p, dt);
    propagate_function(vec_m, dt);

    jacobian.col(i) = (vec_p - vec_m).cast<double>() / (2.0 * eps);
  }
}

static std::ifstream OpenTestDataFile(const std::string& filename) {
  std::filesystem::path lupnt_data_path = GetDataPath();
  std::filesystem::path test_data_path = lupnt_data_path.parent_path() / "test" / "data";
  std::filesystem::path file_path = test_data_path / filename;
  // Open text file
  std::ifstream file(file_path, std::ios::in);
  assert(file.is_open() && "Could not open file");
  return file;
}

static std::vector<double> ReadVector(std::ifstream& file, bool skip_line = false) {
  std::string line;
  if (skip_line) {
    assert(std::getline(file, line));
  }
  assert(std::getline(file, line));
  std::istringstream iss(line);
  std::vector<double> vec;
  double val;
  while (iss >> val) {
    vec.push_back(val);
  }
  return vec;
}

static Vec3 ReadVec3(std::ifstream& file, bool skip_line = false) {
  std::vector<double> vec = ReadVector(file, skip_line);
  return Vec3(vec[0], vec[1], vec[2]);
}
