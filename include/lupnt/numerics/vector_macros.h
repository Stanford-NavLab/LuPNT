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

// Function:
// Real = func(Real)
// New definitions:
// Vec = func(Vec)
#define VEC_DEF_REAL(func) VecX func(const VecX &x);

#define VEC_IMP_REAL(func)               \
  VecX func(const VecX &x) {             \
    VecX out(x.size());                  \
    for (int i = 0; i < x.size(); i++) { \
      out(i) = func(x(i));               \
    }                                    \
    return out;                          \
  }

// Function:
// Vec<size> = func(Vec<size>)
// New definitions:
// Mat<-1,size> = (Mat<-1,size>)
#define VEC_DEF_VECTOR(func, size) Mat<-1, size> func(const Mat<-1, size> &x);

#define VEC_IMP_VECTOR(func, size)             \
  Mat<-1, size> func(const Mat<-1, size> &x) { \
    Mat<-1, size> out(x.rows(), size);         \
    for (int i = 0; i < x.rows(); i++) {       \
      Vec<size> x_ = x.row(i);                 \
      out.row(i) = func(x_);                   \
    }                                          \
    return out;                                \
  }

// Function:
// Vec<size> = func(Vec<size>, Real)
// New definitions:
// Mat<-1,size> = func(Vec<size>, VecX)
// Mat<-1,size> = func(Mat<-1,size>, Real)
// Mat<-1,size> = func(Mat<-1,size>, VecX)
#define VEC_DEF_VECTOR_REAL(func, size)                  \
  Mat<-1, size> func(const Vec<size> &x, const VecX &y); \
  Mat<-1, size> func(const Mat<-1, size> &x, Real y);    \
  Mat<-1, size> func(const Mat<-1, size> &x, const VecX &y);

#define VEC_IMP_VECTOR_REAL(func, size)                         \
  Mat<-1, size> func(const Vec<size> &x, const VecX &y) {       \
    Mat<-1, size> out(y.rows(), size);                          \
    for (int i = 0; i < y.rows(); i++) {                        \
      out.row(i) = func(x, y(i));                               \
    }                                                           \
    return out;                                                 \
  }                                                             \
  Mat<-1, size> func(const Mat<-1, size> &x, Real y) {          \
    Mat<-1, size> out(x.rows(), size);                          \
    for (int i = 0; i < x.rows(); i++) {                        \
      Vec<size> x_ = x.row(i);                                  \
      out.row(i) = func(x_, y);                                 \
    }                                                           \
    return out;                                                 \
  }                                                             \
  Mat<-1, size> func(const Mat<-1, size> &x, const VecX &y) {   \
    ASSERT_WITH_MESSAGE(x.rows() == y.rows(), "Size mismatch"); \
    Mat<-1, size> out(x.rows(), size);                          \
    for (int i = 0; i < x.rows(); i++) {                        \
      Vec<size> x_ = x.row(i);                                  \
      out.row(i) = func(x_, y(i));                              \
    }                                                           \
    return out;                                                 \
  }
// Function:
// Vec<size> = func(Vec<size>, Real, Real)
// New definitions:
// Mat<-1,size> = func(Mat<-1,size>, Real, Real)
#define VEC_DEF_VECTOR_REAL_REAL(func, size) \
  Mat<-1, size> func(const Mat<-1, size> &x, Real y, Real z);

#define VEC_IMP_VECTOR_REAL_REAL(func, size)                   \
  Mat<-1, size> func(const Mat<-1, size> &x, Real y, Real z) { \
    Mat<-1, size> out(x.rows(), size);                         \
    for (int i = 0; i < x.rows(); i++) {                       \
      Vec<size> x_ = x.row(i);                                 \
      out.row(i) = func(x_, y, z);                             \
    }                                                          \
    return out;                                                \
  }

// Function:
// Vec<size> = func(Vec<size>, Vec<size>, Real)
// New definitions:
// Mat<-1,size> = func(Mat<-1,size>, Vec<size>, Real)
// Mat<-1,size> = func(Vec<size>, Mat<-1,size>, Real)
// Mat<-1,size> = func(Mat<-1,size>, Mat<-1,size>, Real)
#define VEC_DEF_VECTOR_VECTOR_REAL(func, size)                            \
  Mat<-1, size> func(const Mat<-1, size> &x, const Vec<size> &y, Real z); \
  Mat<-1, size> func(const Vec<size> &x, const Mat<-1, size> &y, Real z); \
  Mat<-1, size> func(const Mat<-1, size> &x, const Mat<-1, size> &y, Real z);

#define VEC_IMP_VECTOR_VECTOR_REAL(func, size)                                 \
  Mat<-1, size> func(const Mat<-1, size> &x, const Vec<size> &y, Real z) {     \
    Mat<-1, size> out(x.rows(), size);                                         \
    for (int i = 0; i < x.rows(); i++) {                                       \
      Vec<size> x_ = x.row(i);                                                 \
      out.row(i) = func(x_, y, z);                                             \
    }                                                                          \
    return out;                                                                \
  }                                                                            \
  Mat<-1, size> func(const Vec<size> &x, const Mat<-1, size> &y, Real z) {     \
    Mat<-1, size> out(y.rows(), size);                                         \
    for (int i = 0; i < y.rows(); i++) {                                       \
      Vec<size> y_ = y.row(i);                                                 \
      out.row(i) = func(x, y_, z);                                             \
    }                                                                          \
    return out;                                                                \
  }                                                                            \
  Mat<-1, size> func(const Mat<-1, size> &x, const Mat<-1, size> &y, Real z) { \
    ASSERT_WITH_MESSAGE(x.rows() == y.rows(), "Size mismatch");                \
    Mat<-1, size> out(x.rows(), size);                                         \
    for (int i = 0; i < x.rows(); i++) {                                       \
      Vec<size> x_ = x.row(i);                                                 \
      Vec<size> y_ = y.row(i);                                                 \
      out.row(i) = func(x_, y_, z);                                            \
    }                                                                          \
    return out;                                                                \
  }

// Function:
// Vec<size> = func(Vec<size>, Vec<size>
// New definitions:
// Mat<-1,size> = func(Mat<-1,size>, Mat<-1,size>)
// Mat<-1,size> = func(Mat<-1,size>, Vec<size>)
// Mat<-1,size> = func(Vec<size>, Mat<-1,size>)
#define VEC_DEF_VECTOR_VECTOR(func, size)                             \
  Mat<-1, size> func(const Mat<-1, size> &x, const Mat<-1, size> &y); \
  Mat<-1, size> func(const Mat<-1, size> &x, const Vec<size> &y);     \
  Mat<-1, size> func(const Vec<size> &x, const Mat<-1, size> &y);

#define VEC_IMP_VECTOR_VECTOR(func, size)                              \
  Mat<-1, size> func(const Mat<-1, size> &x, const Mat<-1, size> &y) { \
    ASSERT_WITH_MESSAGE(x.rows() == y.rows(), "Size mismatch");        \
    Mat<-1, size> out(x.rows(), size);                                 \
    for (int i = 0; i < x.rows(); i++) {                               \
      Vec<size> x_ = x.row(i);                                         \
      Vec<size> y_ = y.row(i);                                         \
      out.row(i) = func(x_, y_);                                       \
    }                                                                  \
    return out;                                                        \
  }                                                                    \
  Mat<-1, size> func(const Mat<-1, size> &x, const Vec<size> &y) {     \
    Mat<-1, size> out(x.rows(), size);                                 \
    for (int i = 0; i < x.rows(); i++) {                               \
      Vec<size> x_ = x.row(i);                                         \
      out.row(i) = func(x_, y);                                        \
    }                                                                  \
    return out;                                                        \
  }                                                                    \
  Mat<-1, size> func(const Vec<size> &x, const Mat<-1, size> &y) {     \
    Mat<-1, size> out(y.rows(), size);                                 \
    for (int i = 0; i < y.rows(); i++) {                               \
      Vec<size> y_ = y.row(i);                                         \
      out.row(i) = func(x, y_);                                        \
    }                                                                  \
    return out;                                                        \
  }

// Function:
// Real = func(Real, Real)
// New definitions:
// vector = func(vector, Real)
// vector = func(Real, vector)
// vector = func(vector, vector)
#define VEC_DEF_REAL_REAL(func)     \
  VecX func(const VecX &x, Real y); \
  VecX func(Real x, const VecX &y); \
  VecX func(const VecX &x, const VecX &y);

#define VEC_IMP_REAL_REAL(func)                      \
  VecX func(const VecX &x, Real y) {                 \
    VecX out(x.size());                              \
    for (int i = 0; i < x.size(); i++) {             \
      out(i) = func(x(i), y);                        \
    }                                                \
    return out;                                      \
  }                                                  \
  VecX func(Real x, const VecX &y) {                 \
    VecX out(y.size());                              \
    for (int i = 0; i < y.size(); i++) {             \
      out(i) = func(x, y(i));                        \
    }                                                \
    return out;                                      \
  }                                                  \
  VecX func(const VecX &x, const VecX &y) {          \
    assert(x.size() == y.size() && "Size mismatch"); \
    VecX out(x.size());                              \
    for (int i = 0; i < x.size(); i++) {             \
      out(i) = func(x(i), y(i));                     \
    }                                                \
    return out;                                      \
  }
