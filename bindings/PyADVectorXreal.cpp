//                  _  _
//  _   _|_ _  _|o_|__|_
// (_||_||_(_)(_|| |  |
//
// automatic differentiation made easier in C++
// https://github.com/autodiff/autodiff
//
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
//
// Copyright (c) 2018-2022 Allan Leal
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// pybind11 includes
#include "pybind11.hxx"

// autodiff includes
#include <autodiff/forward/real/eigen.hpp>
#include <autodiff/forward/real/real.hpp>

#include "eigen.hxx"
using namespace autodiff;

using Matrix6real = Eigen::Matrix<real, 6, 6, 0, 6, 6>;
using Vector6real = Eigen::Matrix<real, 6, 1, 0, 6, 1>;

void exportVectorXreal0th(py::module& m) {
  exportVector<VectorXreal0th, real0th, isarray(false), isconst(false),
               isview(false)>(m, "VectorXreal0th");
  exportVector<Eigen::Ref<VectorXreal0th>, real0th, isarray(false),
               isconst(false), isview(true)>(m, "VectorXreal0thRef");
  exportVector<Eigen::Ref<const VectorXreal0th>, real0th, isarray(false),
               isconst(true), isview(true)>(m, "VectorXreal0thConstRef");
}

void exportVectorXreal1st(py::module& m) {
  exportVector<VectorXreal1st, real1st, isarray(false), isconst(false),
               isview(false)>(m, "VectorXreal1st");
  exportVector<Eigen::Ref<VectorXreal1st>, real1st, isarray(false),
               isconst(false), isview(true)>(m, "VectorXreal1stRef");
  exportVector<Eigen::Ref<const VectorXreal1st>, real1st, isarray(false),
               isconst(true), isview(true)>(m, "VectorXreal1stConstRef");
}

void exportVectorXreal2nd(py::module& m) {
  exportVector<VectorXreal2nd, real2nd, isarray(false), isconst(false),
               isview(false)>(m, "VectorXreal2nd");
  exportVector<Eigen::Ref<VectorXreal2nd>, real2nd, isarray(false),
               isconst(false), isview(true)>(m, "VectorXreal2ndRef");
  exportVector<Eigen::Ref<const VectorXreal2nd>, real2nd, isarray(false),
               isconst(true), isview(true)>(m, "VectorXreal2ndConstRef");
}

void exportVectorXreal3rd(py::module& m) {
  exportVector<VectorXreal3rd, real3rd, isarray(false), isconst(false),
               isview(false)>(m, "VectorXreal3rd");
  exportVector<Eigen::Ref<VectorXreal3rd>, real3rd, isarray(false),
               isconst(false), isview(true)>(m, "VectorXreal3rdRef");
  exportVector<Eigen::Ref<const VectorXreal3rd>, real3rd, isarray(false),
               isconst(true), isview(true)>(m, "VectorXreal3rdConstRef");
}

void exportVectorXreal4th(py::module& m) {
  exportVector<VectorXreal4th, real4th, isarray(false), isconst(false),
               isview(false)>(m, "VectorXreal4th");
  exportVector<Eigen::Ref<VectorXreal4th>, real4th, isarray(false),
               isconst(false), isview(true)>(m, "VectorXreal4thRef");
  exportVector<Eigen::Ref<const VectorXreal4th>, real4th, isarray(false),
               isconst(true), isview(true)>(m, "VectorXreal4thConstRef");
}

void exportVector2real(py::module& m) {
  exportVector<Vector2real, real, isarray(false), isconst(false),
               isview(false)>(m, "Vector2real");
  exportVector<Eigen::Ref<Vector2real>, real, isarray(false), isconst(false),
               isview(true)>(m, "Vector2realRef");
  exportVector<Eigen::Ref<const Vector2real>, real, isarray(false),
               isconst(true), isview(true)>(m, "Vector2realConstRef");
}

void exportVector3real(py::module& m) {
  exportVector<Vector3real, real, isarray(false), isconst(false),
               isview(false)>(m, "Vector3real");
  exportVector<Eigen::Ref<Vector3real>, real, isarray(false), isconst(false),
               isview(true)>(m, "Vector3realRef");
  exportVector<Eigen::Ref<const Vector3real>, real, isarray(false),
               isconst(true), isview(true)>(m, "Vector3realConstRef");
}

void exportVector4real(py::module& m) {
  exportVector<Vector4real, real, isarray(false), isconst(false),
               isview(false)>(m, "Vector4real");
  exportVector<Eigen::Ref<Vector4real>, real, isarray(false), isconst(false),
               isview(true)>(m, "Vector4realRef");
  exportVector<Eigen::Ref<const Vector4real>, real, isarray(false),
               isconst(true), isview(true)>(m, "Vector4realConstRef");
}

void exportVector6real(py::module& m) {
  exportVector<Vector6real, real, isarray(false), isconst(false),
               isview(false)>(m, "Vector6real");
  exportVector<Eigen::Ref<Vector6real>, real, isarray(false), isconst(false),
               isview(true)>(m, "Vector6realRef");
  exportVector<Eigen::Ref<const Vector6real>, real, isarray(false),
               isconst(true), isview(true)>(m, "Vector6realConstRef");
}
