#pragma once

#include <matplot/matplot.h>
#include <omp.h>
#include <spdlog/spdlog.h>
#include <yaml-cpp/yaml.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <highfive/H5Easy.hpp>
#include <iostream>
#include <magic_enum/magic_enum.hpp>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unsupported/Eigen/EulerAngles>
#include <unsupported/Eigen/SpecialFunctions>
#include <vector>

#define DEFINE_VECTORS_MATRICES_FIXED(size)                   \
  using Vec##size = Matrix<Real, size, 1>;                    \
  using Vec##size##d = Matrix<double, size, 1>;               \
  using Vec##size##i = Matrix<int, size, 1>;                  \
  using Vec##size##b = Matrix<bool, size, 1>;                 \
  using Mat##size = Matrix<Real, size, size>;                 \
  using Mat##size##d = Matrix<double, size, size>;            \
  using Mat##size##i = Matrix<int, size, size>;               \
  using Mat##size##b = Matrix<bool, size, size>;              \
  using RowVec##size = Matrix<Real, 1, size>;                 \
  using RowVec##size##d = Matrix<double, 1, size>;            \
  using RowVec##size##i = Matrix<int, 1, size>;               \
  using RowVec##size##b = Matrix<bool, 1, size>;              \
  using MatX##size = Matrix<Real, Eigen::Dynamic, size>;      \
  using MatX##size##d = Matrix<double, Eigen::Dynamic, size>; \
  using MatX##size##i = Matrix<int, Eigen::Dynamic, size>;    \
  using MatX##size##b = Matrix<bool, Eigen::Dynamic, size>;   \
  using Mat##size##X = Matrix<Real, size, Eigen::Dynamic>;    \
  using Mat##size##Xd = Matrix<double, size, Eigen::Dynamic>; \
  using Mat##size##Xi = Matrix<int, size, Eigen::Dynamic>;    \
  using Mat##size##Xb = Matrix<bool, size, Eigen::Dynamic>;   \
  using Arr##size = Array<Real, size, 1>;                     \
  using Arr##size##d = Array<double, size, 1>;                \
  using Arr##size##i = Array<int, size, 1>;                   \
  using Arr##size##b = Array<bool, size, 1>;                  \
  using Arr##size##size = Array<Real, size, size>;            \
  using Arr##size##size##d = Array<double, size, size>;       \
  using Arr##size##size##i = Array<int, size, size>;          \
  using Arr##size##size##b = Array<bool, size, size>;         \
  using ArrX##size = Array<Real, Eigen::Dynamic, size>;       \
  using ArrX##size##d = Array<double, Eigen::Dynamic, size>;  \
  using ArrX##size##i = Array<int, Eigen::Dynamic, size>;     \
  using ArrX##size##b = Array<bool, Eigen::Dynamic, size>;    \
  using Arr##size##X = Array<Real, size, Eigen::Dynamic>;     \
  using Arr##size##Xd = Array<double, size, Eigen::Dynamic>;  \
  using Arr##size##Xi = Array<int, size, Eigen::Dynamic>;     \
  using Arr##size##Xb = Array<bool, size, Eigen::Dynamic>;

#define DEFINE_VECTORS_MATRICES_DYNAMIC()                       \
  using VecX = Matrix<Real, Eigen::Dynamic, 1>;                 \
  using VecXd = Matrix<double, Eigen::Dynamic, 1>;              \
  using VecXi = Matrix<int, Eigen::Dynamic, 1>;                 \
  using VecXb = Matrix<bool, Eigen::Dynamic, 1>;                \
  using MatX = Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;    \
  using MatXd = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; \
  using MatXi = Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;    \
  using MatXb = Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;   \
  using RowVecX = Matrix<Real, 1, Eigen::Dynamic>;              \
  using RowVecXd = Matrix<double, 1, Eigen::Dynamic>;           \
  using RowVecXi = Matrix<int, 1, Eigen::Dynamic>;              \
  using RowVecXb = Matrix<bool, 1, Eigen::Dynamic>;             \
  using ArrX = Array<Real, Eigen::Dynamic, Eigen::Dynamic>;     \
  using ArrXd = Array<double, Eigen::Dynamic, Eigen::Dynamic>;  \
  using ArrXi = Array<int, Eigen::Dynamic, Eigen::Dynamic>;     \
  using ArrXb = Array<bool, Eigen::Dynamic, Eigen::Dynamic>;

#define DEFINE_VECTORS_MATRICES()   \
  DEFINE_VECTORS_MATRICES_FIXED(1)  \
  DEFINE_VECTORS_MATRICES_FIXED(2)  \
  DEFINE_VECTORS_MATRICES_FIXED(3)  \
  DEFINE_VECTORS_MATRICES_FIXED(4)  \
  DEFINE_VECTORS_MATRICES_FIXED(5)  \
  DEFINE_VECTORS_MATRICES_FIXED(6)  \
  DEFINE_VECTORS_MATRICES_FIXED(7)  \
  DEFINE_VECTORS_MATRICES_FIXED(8)  \
  DEFINE_VECTORS_MATRICES_FIXED(9)  \
  DEFINE_VECTORS_MATRICES_FIXED(10) \
  DEFINE_VECTORS_MATRICES_DYNAMIC()

namespace lupnt {

  // Pointers
  template <typename T> using Ptr = std::shared_ptr<T>;
  template <typename T> using UniquePtr = std::unique_ptr<T>;
  template <typename T> using SharedPtr = std::shared_ptr<T>;

  /// @brief Construct a `Ptr<T>` (alias for `std::shared_ptr<T>`), forwarding
  /// constructor arguments to `T`. Used throughout LuPNT to create
  /// Agent/Device/Dynamics/Filter/etc. instances managed by shared ownership.
  template <typename T, typename... Args> Ptr<T> MakePtr(Args&&... args) {
    return std::make_shared<T>(std::forward<Args>(args)...);
  }

  /// @brief Construct a `UniquePtr<T>` (alias for `std::unique_ptr<T>`), forwarding
  /// constructor arguments to `T`.
  template <typename T, typename... Args> UniquePtr<T> MakeUnique(Args&&... args) {
    return std::make_unique<T>(std::forward<Args>(args)...);
  }

  /// @brief Construct a `SharedPtr<T>` (alias for `std::shared_ptr<T>`), forwarding
  /// constructor arguments to `T`.
  template <typename T, typename... Args> SharedPtr<T> MakeShared(Args&&... args) {
    return std::make_shared<T>(std::forward<Args>(args)...);
  }

  // Eigen
  using Eigen::Array;
  using Eigen::Block;
  using Eigen::DenseBase;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixBase;
  using Eigen::MatrixX;
  using Eigen::Vector;
  using Eigen::VectorX;

  namespace ad = autodiff;

  using Real = ad::real;
  using ad::at;
  using ad::derivative;
  using ad::jacobian;
  using ad::wrt;
  template <int rows, int cols> using Mat = Matrix<Real, rows, cols>;
  template <int rows, int cols> using Matd = Matrix<double, rows, cols>;
  template <typename T, int size> using VecT = Matrix<T, size, 1>;
  template <int size> using Vec = Matrix<Real, size, 1>;
  template <int size> using Vecd = Matrix<double, size, 1>;
  template <int size> using RowVecd = Matrix<double, 1, size>;
  using RowVecXd = Matrix<double, 1, Eigen::Dynamic>;

  DEFINE_VECTORS_MATRICES()

  // Print Formats
  static Eigen::IOFormat FMT_CLEAN(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
  static Eigen::IOFormat FMT_HEAVY(Eigen::FullPrecision, 0, ", ", ",\n", "[", "]", "[", "]");
  static Eigen::IOFormat FMT_COMPACT(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ";\n ", "",
                                     "", "[", "]");
  // Magic Enum
  using magic_enum::enum_cast;
  using magic_enum::enum_integer;
  using magic_enum::enum_name;
  using magic_enum::enum_names;
  using magic_enum::enum_value;
  using magic_enum::enum_values;

  using json = nlohmann::json;
  using Config = YAML::Node;

  /// @brief Stream a `std::vector<T>` as a bracketed, comma-separated list (e.g.
  /// `[1, 2, 3]`). Used by `Logger`/`fmt`-based debug output throughout LuPNT to
  /// print vectors of bodies, names, etc. in log messages.
  template <typename T> std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os << "[";
    for (size_t i = 0; i < vec.size(); i++) {
      os << vec[i];
      if (i != vec.size() - 1) os << ", ";
    }
    os << "]";
    return os;
  }

  /// @brief Get the global simulation reference epoch [s, TDB since J2000].
  ///
  /// This is the absolute time corresponding to simulation time `t = 0`. It is set once
  /// by `Simulation`'s constructor from the config's `epoch` field (via `SetLupntEpoch`)
  /// and then read throughout the dynamics/agents/conversions modules (e.g.
  /// `NumericalOrbitDynamics::CalcContrib`, `AttitudeDynamics`, `Satellite`,
  /// `Constellation`/TLE propagation, `CesiumViewer`) to convert a local simulation time
  /// `t` into an absolute TDB epoch via `t_tdb = GetLupntEpoch() + t`.
  ///
  /// @return Reference epoch [s, TDB since J2000]
  Real GetLupntEpoch();

  /// @brief Set the global simulation reference epoch [s, TDB since J2000].
  ///
  /// Called once during `Simulation` construction (from the config's `epoch` field) to
  /// anchor simulation time `t = 0` to an absolute TDB epoch; also used directly in
  /// example scenarios that set up a propagation without a full `Simulation` config.
  ///
  /// @param epoch Reference epoch to install [s, TDB since J2000]
  void SetLupntEpoch(Real epoch);

}  // namespace lupnt

template <> struct fmt::formatter<lupnt::Real> : fmt::formatter<double> {
  template <typename FormatContext> auto format(const lupnt::Real& x, FormatContext& ctx) const {
    return fmt::formatter<double>::format(static_cast<double>(x), ctx);
  }
};

// Formatter for all Eigen types (matrices, vectors, expressions like Transpose, etc.)
template <typename Derived> struct fmt::formatter<Eigen::EigenBase<Derived>>
    : fmt::formatter<std::string> {
  template <typename FormatContext>
  auto format(const Eigen::EigenBase<Derived>& mat, FormatContext& ctx) const {
    std::ostringstream oss;
    // Check if it's a DenseBase type (has .format() method)
    if constexpr (std::is_base_of_v<Eigen::DenseBase<Derived>, Derived>) {
      oss << static_cast<const Eigen::DenseBase<Derived>&>(mat).format(lupnt::FMT_COMPACT);
    } else {
      oss << mat.derived();
    }
    return fmt::formatter<std::string>::format(oss.str(), ctx);
  }
};

// Specific formatter for Eigen::Transpose
template <typename MatrixType> struct fmt::formatter<Eigen::Transpose<MatrixType>>
    : fmt::formatter<std::string> {
  template <typename FormatContext>
  auto format(const Eigen::Transpose<MatrixType>& mat, FormatContext& ctx) const {
    std::ostringstream oss;
    oss << mat.eval().format(lupnt::FMT_COMPACT);
    return fmt::formatter<std::string>::format(oss.str(), ctx);
  }
};
