// Microbenchmark comparing Jacobian computation through autodiff, Eigen's
// AutoDiffJacobian, and numerical differentiation.
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/NumericalDiff>
#include <vector>

namespace ad = autodiff;
using ad::at;
using ad::dual;
using ad::jacobian;
using ad::real;
using ad::wrt;
using Eigen::Matrix;
using Eigen::Vector;
using namespace std::chrono;

constexpr int N_FUNC = 1'000;
constexpr int N_CALLS = 100'000;

template <typename T, int N> Vector<T, N> MyFunction(const Vector<T, N>& x) {
  // Deliberately expensive nonlinear function so call overhead does not dominate
  // every benchmark case.
  Vector<T, N> y;
  if constexpr (N == Eigen::Dynamic) y.resize(2);
  y.setZero();

  for (int i = 0; i < N_FUNC; ++i) {
    for (int j = 0; j < y.size(); ++j) {
      y(j) += sin(x(j) + i) * exp(x((j + 1) % y.size()) - i)
              + log(1 + x(j) * x((j + 1) % y.size()) + i)
              + cos(x(j) - i) * (x((j + 1) % y.size()) * x((j + 1) % y.size()) + i);
    }
  }
  return y;
}

template <typename T, int N> void AutodiffJacobian(int num_calls) {
  auto func = [](const Vector<T, N>& x) { return MyFunction<T, N>(x); };

  Vector<T, N> x, y;
  Matrix<double, N, N> J;
  if constexpr (N == Eigen::Dynamic) {
    x.resize(2);
    y.resize(2);
    J.resize(2, 2);
  }
  x << 1.0, 2.0;

  auto start = high_resolution_clock::now();
  for (int i = 0; i < num_calls; ++i) {
    jacobian(func, wrt(x), at(x), y, J);
  }
  auto end = high_resolution_clock::now();

  std::cout << "y = " << y.transpose() << "\n";
  std::cout << "J = \n" << J << "\n";
  std::cout << std::fixed << std::setprecision(3)
            << duration<double, std::milli>(end - start).count() / 1000 << " s\n\n";
}

template <int N, bool Numerical = false> struct EigenFunctor {
  using Scalar = double;
  using InputType = Matrix<double, N, 1>;
  using ValueType = Matrix<double, N, 1>;
  using JacobianType = Matrix<double, N, N>;
  static constexpr int InputsAtCompileTime = N;
  static constexpr int ValuesAtCompileTime = N;

  template <typename ScalarT>
  void operator()(const Matrix<ScalarT, N, 1>& x, Matrix<ScalarT, N, 1>* y) const {
    *y = MyFunction<ScalarT, N>(x);
  }
  void operator()(const InputType& x, ValueType& y) const { y = MyFunction<double, N>(x); }
  constexpr int values() const { return 2; }
};

template <int N> void EigenAutoDiffJacobian(int num_calls) {
  EigenFunctor<N, false> func;
  Eigen::AutoDiffJacobian<EigenFunctor<N, false>> autodiffFunc(func);

  Vector<double, N> x, y;
  Matrix<double, N, N> J;
  if constexpr (N == Eigen::Dynamic) {
    x.resize(2);
    y.resize(2);
    J.resize(2, 2);
  }
  x << 1.0, 2.0;

  auto start = high_resolution_clock::now();
  for (int i = 0; i < num_calls; ++i) {
    autodiffFunc(x, &y, &J);
  }
  auto end = high_resolution_clock::now();

  std::cout << "y = " << y.transpose() << "\n";
  std::cout << "J = \n"
            << J << "\n"
            << std::fixed << std::setprecision(3)
            << "t = " << duration<double, std::milli>(end - start).count() / 1000 << " s\n\n";
}

template <int N> void EigenNumericalJacobian(int num_calls) {
  EigenFunctor<N, true> func;
  Eigen::NumericalDiff<EigenFunctor<N, true>, Eigen::Central> numDiff(func);

  Vector<double, N> x, y;
  Matrix<double, N, N> J;
  if constexpr (N == Eigen::Dynamic) {
    x.resize(2);
    y.resize(2);
    J.resize(2, 2);
  }
  x << 1.0, 2.0;

  auto start = high_resolution_clock::now();
  for (int i = 0; i < num_calls; ++i) {
    func(x, y);
    numDiff.df(x, J);
  }
  auto end = high_resolution_clock::now();

  std::cout << "y = " << y.transpose() << "\n";
  std::cout << "J = \n"
            << J << "\n"
            << std::fixed << std::setprecision(3)
            << "t = " << duration<double, std::milli>(end - start).count() / 1000 << " s\n\n";
}

template <int N> void FunctionCalls(int num_calls) {
  Matrix<double, N, 1> x, y;
  if constexpr (N == Eigen::Dynamic) x.resize(2);
  x << 1.0, 2.0;

  auto start = high_resolution_clock::now();
  for (int i = 0; i < num_calls; ++i) {
    y = MyFunction<double, N>(x);
  }
  auto end = high_resolution_clock::now();

  std::cout << "y = " << y.transpose() << "\n";
  std::cout << std::fixed << std::setprecision(3)
            << "t = " << duration<double, std::milli>(end - start).count() / 1000 << " s\n\n";
}

int main() {
  std::cout << "Double fixed" << std::endl;
  FunctionCalls<2>(N_CALLS);

  std::cout << "Double dynamic" << std::endl;
  FunctionCalls<Eigen::Dynamic>(N_CALLS);

  std::cout << "Dual fixed" << std::endl;
  AutodiffJacobian<dual, 2>(N_CALLS);

  std::cout << "Dual dynamic" << std::endl;
  AutodiffJacobian<dual, Eigen::Dynamic>(N_CALLS);

  std::cout << "Real fixed" << std::endl;
  AutodiffJacobian<real, 2>(N_CALLS);

  std::cout << "Real dynamic" << std::endl;
  AutodiffJacobian<real, Eigen::Dynamic>(N_CALLS);

  std::cout << "Eigen autodiff fixed" << std::endl;
  EigenAutoDiffJacobian<2>(N_CALLS);

  std::cout << "Eigen autodiff dynamic" << std::endl;
  EigenAutoDiffJacobian<Eigen::Dynamic>(N_CALLS);

  std::cout << "Eigen numerical fixed" << std::endl;
  EigenNumericalJacobian<2>(N_CALLS);

  std::cout << "Eigen numerical dynamic" << std::endl;
  EigenNumericalJacobian<Eigen::Dynamic>(N_CALLS);

  return 0;
}
