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
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
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

// matrix utils
#include "Matrix.hpp"

// C++ includes
#include <sstream>
#include <vector>

// autodiff includes
#include <autodiff/common/meta.hpp>
#include <autodiff/forward/real/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <autodiff/forward/utils/gradient.hpp>
using namespace autodiff;
using autodiff::detail::isSame;

template<typename Fun, typename T>
void exportGradient(py::module& m, const char* typestr)
{
    m.def(typestr, [](const Fun& fun, T x)
    {
        return gradient(fun, wrt(x), at(x));
        
    }, py::arg("func"), py::arg("wrt"), py::return_value_policy::reference);
}

template<typename Fun, typename T>
void exportJacobian(py::module& m, const char* typestr)
{
    m.def(typestr, [](const Fun& fun, T x)
    {
        const Eigen::MatrixXd J = jacobian(fun, wrt(x), at(x));
        return EigenMatrixToPythonArray(J);   // we cannot directly use MatriXd since we cannot import <pybind/eigen.h> in this project -> convert to Python array
        
    }, py::arg("func"), py::arg("wrt"), py::return_value_policy::reference);
}

void export_gradient(py::module& m) { exportGradient<std::function<real(VectorXreal)>, VectorXreal>(m, "gradient"); }
void export_jacobian(py::module& m) { exportJacobian<std::function<VectorXreal(VectorXreal)>, VectorXreal>(m, "jacobian"); }