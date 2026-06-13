/**
 * @file orbit_state_utils.h
 * @author Stanford NAV LAB
 * @brief Util functions for state conversions
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <functional>
#include <map>

#include "lupnt/core/constants.h"

/// @brief Declare/define a `VecX func(const VecX& x)` overload that applies a scalar
/// `Real func(Real)` element-wise.
///
/// `VEC_DEF_REAL` declares the overload (used in headers, e.g. time_conversions.h for
/// `UtcToUt1`, `TaiToTt`, etc.); `VEC_IMP_REAL` defines it (used in the corresponding .cc).
/// Lets time/anomaly conversion functions be applied to a whole vector of epochs/angles at
/// once (in an OpenMP-parallel loop), e.g. converting an entire epoch history from UTC to
/// TAI.
// Function:
// Real = func(Real)
// New definitions:
// Vec = func(Vec)
#define VEC_DEF_REAL(func) VecX func(const VecX& x);

#define VEC_IMP_REAL(func)                                                                  \
  VecX func(const VecX& x) {                                                                \
    VecX out(x.size());                                                                     \
    _Pragma("omp parallel for") for (int i = 0; i < x.size(); i++) { out(i) = func(x(i)); } \
    return out;                                                                             \
  }

/// @brief Declare/define a `Vec<size> func(const Vec<size>&)` overload (delegating to the
/// `State`-based implementation) plus a `Mat<-1,size> func(const Mat<-1,size>&)` overload
/// that applies it row-wise.
///
/// `VEC_DEF_VECTOR` declares both overloads (e.g. attitude_conversions.h's
/// `ScalarFirstToLast`/`ScalarLastToFirst`/`NormalizeQuat` for quaternions,
/// coordinate_conversions.h's `EastNorthUpToAzElRange`/`AzElRangeToEastNorthUp`);
/// `VEC_IMP_VECTOR` defines them. The matrix overload lets a fixed-size state/vector
/// conversion be applied to a batch of states (one per row), e.g. converting a whole
/// trajectory's quaternions in one call.
// Function:
// Vec<size> = func(Vec<size>)
// New definitions:
// Mat<-1,size> = (Mat<-1,size>)
#define VEC_DEF_VECTOR(func, size)    \
  Vec<size> func(const Vec<size>& x); \
  Mat<-1, size> func(const Mat<-1, size>& x);

#define VEC_IMP_VECTOR(func, size)                                                  \
  Vec<size> func(const Vec<size>& x) { return func(static_cast<const State&>(x)); } \
  Mat<-1, size> func(const Mat<-1, size>& x) {                                      \
    Mat<-1, size> out(x.rows(), size);                                              \
    _Pragma("omp parallel for") for (int i = 0; i < x.rows(); i++) {                \
      Vec<size> x_ = x.row(i);                                                      \
      out.row(i) = func(x_);                                                        \
    }                                                                               \
    return out;                                                                     \
  }

/// @brief Like `VEC_DEF_VECTOR`/`VEC_IMP_VECTOR`, but for conversions that change the
/// state's dimensionality: `Vec<size1> -> Vec<size2>` (and the corresponding
/// `Mat<-1,size1> -> Mat<-1,size2>` row-wise batch overload).
// Function:
// Vec<size> = func(Vec<size>)
// New definitions:
// Mat<-1,size> = (Mat<-1,size>)
#define VEC_DEF_VECTOR_SIZE(func, size1, size2) \
  Vec<size2> func(const Vec<size1>& x);         \
  Mat<-1, size2> func(const Mat<-1, size1>& x);

#define VEC_IMP_VECTOR_SIZE(func, size1, size2)                                       \
  Vec<size2> func(const Vec<size1>& x) { return func(static_cast<const State&>(x)); } \
  Mat<-1, size2> func(const Mat<-1, size1>& x) {                                      \
    Mat<-1, size2> out(x.rows(), size2);                                              \
    _Pragma("omp parallel for") for (int i = 0; i < x.rows(); i++) {                  \
      Vec<size1> x_ = x.row(i);                                                       \
      out.row(i) = func(x_);                                                          \
    }                                                                                 \
    return out;                                                                       \
  }

/// @brief Declare/define batch overloads of a `Vec<size> func(const Vec<size>&, Real)`
/// conversion (e.g. classical orbital elements <-> Cartesian state, which additionally
/// depend on a scalar parameter such as GM): a fixed state with a vector of scalars
/// (`Mat<-1,size> func(Vec<size>, VecX)`), a batch of states with one shared scalar
/// (`Mat<-1,size> func(Mat<-1,size>, Real)`), and a batch of states each paired with its
/// own scalar (`Mat<-1,size> func(Mat<-1,size>, VecX)`).
///
/// `VEC_DEF_VECTOR_REAL` declares these (e.g. state_conversions.h's `ClassicalToCart`,
/// `CartToClassical`, `ClassicalToEquinoctial`, ...); `VEC_IMP_VECTOR_REAL` defines them.
// Function:
// Vec<size> = func(Vec<size>, Real)
// New definitions:
// Mat<-1,size> = func(Vec<size>, VecX)
// Mat<-1,size> = func(Mat<-1,size>, Real)
// Mat<-1,size> = func(Mat<-1,size>, VecX)
#define VEC_DEF_VECTOR_REAL(func, size)                  \
  Vec<size> func(const Vec<size>& x, Real y);            \
  Mat<-1, size> func(const Vec<size>& x, const VecX& y); \
  Mat<-1, size> func(const Mat<-1, size>& x, Real y);    \
  Mat<-1, size> func(const Mat<-1, size>& x, const VecX& y);

#define VEC_IMP_VECTOR_REAL(func, size)                                                        \
  Vec<size> func(const Vec<size>& x, Real y) { return func(static_cast<const State&>(x), y); } \
  Mat<-1, size> func(const Vec<size>& x, const VecX& y) {                                      \
    Mat<-1, size> out(y.rows(), size);                                                         \
    _Pragma("omp parallel for") for (int i = 0; i < y.rows(); i++) out.row(i) = func(x, y(i)); \
    return out;                                                                                \
  }                                                                                            \
  Mat<-1, size> func(const Mat<-1, size>& x, Real y) {                                         \
    Mat<-1, size> out(x.rows(), size);                                                         \
    _Pragma("omp parallel for") for (int i = 0; i < x.rows(); i++) {                           \
      Vec<size> x_ = x.row(i);                                                                 \
      out.row(i) = func(x_, y);                                                                \
    }                                                                                          \
    return out;                                                                                \
  }                                                                                            \
  Mat<-1, size> func(const Mat<-1, size>& x, const VecX& y) {                                  \
    LUPNT_CHECK(x.rows() == y.rows(), "Size mismatch", "VectorMacros");                        \
    Mat<-1, size> out(x.rows(), size);                                                         \
    _Pragma("omp parallel for") for (int i = 0; i < x.rows(); i++) {                           \
      Vec<size> x_ = x.row(i);                                                                 \
      out.row(i) = func(x_, y(i));                                                             \
    }                                                                                          \
    return out;                                                                                \
  }

/// @brief Declare/define a `Mat<-1,size> func(const Mat<-1,size>&, Real, Real)` batch
/// overload of `Vec<size> func(const Vec<size>&, Real, Real)`, applying the conversion
/// row-wise with the same two scalar parameters shared across the whole batch.
///
/// `VEC_DEF_VECTOR_REAL_REAL` declares it, `VEC_IMP_VECTOR_REAL_REAL` defines it.
// Function:
// Vec<size> = func(Vec<size>, Real, Real)
// New definitions:
// Mat<-1,size> = func(Mat<-1,size>, Real, Real)
#define VEC_DEF_VECTOR_REAL_REAL(func, size)          \
  Vec<size> func(const Vec<size>& x, Real y, Real z); \
  Mat<-1, size> func(const Mat<-1, size>& x, Real y, Real z);

#define VEC_IMP_VECTOR_REAL_REAL(func, size)                         \
  Vec<size> func(const Vec<size>& x, Real y, Real z) {               \
    return func(static_cast<const State&>(x), y, z);                 \
  }                                                                  \
  Mat<-1, size> func(const Mat<-1, size>& x, Real y, Real z) {       \
    Mat<-1, size> out(x.rows(), size);                               \
    _Pragma("omp parallel for") for (int i = 0; i < x.rows(); i++) { \
      Vec<size> x_ = x.row(i);                                       \
      out.row(i) = func(x_, y, z);                                   \
    }                                                                \
    return out;                                                      \
  }

/// @brief Declare/define batch overloads of a `Vec<size> func(const Vec<size>&, const
/// Vec<size>&, Real)` conversion (two state-sized vector inputs plus a scalar parameter):
/// a batch of first arguments against a fixed second argument, a fixed first argument
/// against a batch of second arguments, and batch-vs-batch (row-wise).
///
/// `VEC_DEF_VECTOR_VECTOR_REAL` declares these, `VEC_IMP_VECTOR_VECTOR_REAL` defines them.
// Function:
// Vec<size> = func(Vec<size>, Vec<size>, Real)
// New definitions:
// Mat<-1,size> = func(Mat<-1,size>, Vec<size>, Real)
// Mat<-1,size> = func(Vec<size>, Mat<-1,size>, Real)
// Mat<-1,size> = func(Mat<-1,size>, Mat<-1,size>, Real)
#define VEC_DEF_VECTOR_VECTOR_REAL(func, size)                            \
  Mat<-1, size> func(const Mat<-1, size>& x, const Vec<size>& y, Real z); \
  Mat<-1, size> func(const Vec<size>& x, const Mat<-1, size>& y, Real z); \
  Mat<-1, size> func(const Mat<-1, size>& x, const Mat<-1, size>& y, Real z);

#define VEC_IMP_VECTOR_VECTOR_REAL(func, size)                                 \
  Mat<-1, size> func(const Mat<-1, size>& x, const Vec<size>& y, Real z) {     \
    Mat<-1, size> out(x.rows(), size);                                         \
    _Pragma("omp parallel for") for (int i = 0; i < x.rows(); i++) {           \
      Vec<size> x_ = x.row(i);                                                 \
      out.row(i) = func(x_, y, z);                                             \
    }                                                                          \
    return out;                                                                \
  }                                                                            \
  Mat<-1, size> func(const Vec<size>& x, const Mat<-1, size>& y, Real z) {     \
    Mat<-1, size> out(y.rows(), size);                                         \
    _Pragma("omp parallel for") for (int i = 0; i < y.rows(); i++) {           \
      Vec<size> y_ = y.row(i);                                                 \
      out.row(i) = func(x, y_, z);                                             \
    }                                                                          \
    return out;                                                                \
  }                                                                            \
  Mat<-1, size> func(const Mat<-1, size>& x, const Mat<-1, size>& y, Real z) { \
    LUPNT_CHECK(x.rows() == y.rows(), "Size mismatch", "VectorMacros");        \
    Mat<-1, size> out(x.rows(), size);                                         \
    _Pragma("omp parallel for") for (int i = 0; i < x.rows(); i++) {           \
      Vec<size> x_ = x.row(i);                                                 \
      Vec<size> y_ = y.row(i);                                                 \
      out.row(i) = func(x_, y_, z);                                            \
    }                                                                          \
    return out;                                                                \
  }

/// @brief Declare/define batch overloads of a `Vec<size> func(const Vec<size>&, const
/// Vec<size>&)` conversion (two state-sized vector inputs, e.g.
/// `InertialToSynodic`/`SynodicToInertial` in state_conversions.h): batch-vs-batch
/// (row-wise), batch-vs-fixed, and fixed-vs-batch.
///
/// `VEC_DEF_VECTOR_VECTOR` declares these, `VEC_IMP_VECTOR_VECTOR` defines them.
// Function:
// Vec<size> = func(Vec<size>, Vec<size>
// New definitions:
// Mat<-1,size> = func(Mat<-1,size>, Mat<-1,size>)
// Mat<-1,size> = func(Mat<-1,size>, Vec<size>)
// Mat<-1,size> = func(Vec<size>, Mat<-1,size>)
#define VEC_DEF_VECTOR_VECTOR(func, size)                             \
  Vec<size> func(const Vec<size>& x, const Vec<size>& y);             \
  Mat<-1, size> func(const Mat<-1, size>& x, const Mat<-1, size>& y); \
  Mat<-1, size> func(const Mat<-1, size>& x, const Vec<size>& y);     \
  Mat<-1, size> func(const Vec<size>& x, const Mat<-1, size>& y);

#define VEC_IMP_VECTOR_VECTOR(func, size)                                    \
  Vec<size> func(const Vec<size>& x, const Vec<size>& y) {                   \
    return func(static_cast<const State&>(x), static_cast<const State&>(y)); \
  }                                                                          \
  Mat<-1, size> func(const Mat<-1, size>& x, const Mat<-1, size>& y) {       \
    LUPNT_CHECK(x.rows() == y.rows(), "Size mismatch", "VectorMacros");      \
    Mat<-1, size> out(x.rows(), size);                                       \
    _Pragma("omp parallel for") for (int i = 0; i < x.rows(); i++) {         \
      Vec<size> x_ = x.row(i);                                               \
      Vec<size> y_ = y.row(i);                                               \
      out.row(i) = func(x_, y_);                                             \
    }                                                                        \
    return out;                                                              \
  }                                                                          \
  Mat<-1, size> func(const Mat<-1, size>& x, const Vec<size>& y) {           \
    Mat<-1, size> out(x.rows(), size);                                       \
    _Pragma("omp parallel for") for (int i = 0; i < x.rows(); i++) {         \
      Vec<size> x_ = x.row(i);                                               \
      out.row(i) = func(x_, y);                                              \
    }                                                                        \
    return out;                                                              \
  }                                                                          \
  Mat<-1, size> func(const Vec<size>& x, const Mat<-1, size>& y) {           \
    Mat<-1, size> out(y.rows(), size);                                       \
    _Pragma("omp parallel for") for (int i = 0; i < y.rows(); i++) {         \
      Vec<size> y_ = y.row(i);                                               \
      out.row(i) = func(x, y_);                                              \
    }                                                                        \
    return out;                                                              \
  }

/// @brief Declare/define batch overloads of a scalar `Real func(Real, Real)` conversion
/// (e.g. anomaly_conversions.h's `EccToTrueAnomaly`, `MeanToEccAnomaly`, and
/// `GetOrbitalPeriod`, which take an anomaly/semi-major-axis plus eccentricity/GM):
/// vector-vs-scalar, scalar-vs-vector, and element-wise vector-vs-vector.
///
/// `VEC_DEF_REAL_REAL` declares these (throws on size mismatch in the vector-vs-vector
/// case), `VEC_IMP_REAL_REAL` defines them.
// Function:
// Real = func(Real, Real)
// New definitions:
// vector = func(vector, Real)
// vector = func(Real, vector)
// vector = func(vector, vector)
#define VEC_DEF_REAL_REAL(func)     \
  VecX func(const VecX& x, Real y); \
  VecX func(Real x, const VecX& y); \
  VecX func(const VecX& x, const VecX& y);

#define VEC_IMP_REAL_REAL(func)                                                                   \
  VecX func(const VecX& x, Real y) {                                                              \
    VecX out(x.size());                                                                           \
    _Pragma("omp parallel for") for (int i = 0; i < x.size(); i++) { out(i) = func(x(i), y); }    \
    return out;                                                                                   \
  }                                                                                               \
  VecX func(Real x, const VecX& y) {                                                              \
    VecX out(y.size());                                                                           \
    _Pragma("omp parallel for") for (int i = 0; i < y.size(); i++) { out(i) = func(x, y(i)); }    \
    return out;                                                                                   \
  }                                                                                               \
  VecX func(const VecX& x, const VecX& y) {                                                       \
    if (x.size() != y.size()) throw std::runtime_error("Size mismatch");                          \
    VecX out(x.size());                                                                           \
    _Pragma("omp parallel for") for (int i = 0; i < x.size(); i++) { out(i) = func(x(i), y(i)); } \
    return out;                                                                                   \
  }
