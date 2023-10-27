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
#include "../pybind11.hxx"

// C++ includes
#include <sstream>
#include <vector>

// autodiff includes
#include <autodiff/common/meta.hpp>
#include <autodiff/forward/real/real.hpp>
#include <autodiff/forward/utils/derivative.hpp>
using namespace autodiff;
using autodiff::detail::isSame;


template<typename Fun, typename T>
void exportDerivative1(py::module& m, const char* typestr)
{
    m.def(typestr, [](const Fun& fun, T x)
    {
        return derivative(fun, wrt(x), at(x));
    }, py::return_value_policy::reference);
}

template<typename Fun, typename T>
void exportDerivative2(py::module& m, const char* typestr)
{
    m.def(typestr, [](const Fun& fun, T x, T y, int wrtidx)
    {
        switch (wrtidx)
        {
            case 0: return derivative(fun, wrt(x), at(x, y));
            case 1: return derivative(fun, wrt(y), at(x, y));
            default: throw std::invalid_argument("wrtidx must be 0 or 1");
        }
    }, py::arg("func"), py::arg("x"), py::arg("y"), py::arg("wrt"), py::return_value_policy::reference);
}

template<typename Fun, typename T>
void exportDerivative3(py::module& m, const char* typestr)
{
    m.def(typestr, [](const Fun& fun, T x, T y, T z, int wrtidx)
    {
        switch (wrtidx)
        {
            case 0: return derivative(fun, wrt(x), at(x, y, z));
            case 1: return derivative(fun, wrt(y), at(x, y, z));
            case 2: return derivative(fun, wrt(z), at(x, y, z));
            default: throw std::invalid_argument("wrtidx must be 0, 1 or 2");
        }
    }, py::arg("func"), py::arg("x"), py::arg("y"), py::arg("z"), py::arg("wrt"), py::return_value_policy::reference);
}

template<typename Fun, typename T>
void exportDerivative4(py::module& m, const char* typestr)
{
    m.def(typestr, [](const Fun& fun, T w, T x, T y, T z, int wrtidx)
    {
        switch (wrtidx)
        {
            case 0: return derivative(fun, wrt(w), at(w, x, y, z));
            case 1: return derivative(fun, wrt(x), at(w, x, y, z));
            case 2: return derivative(fun, wrt(y), at(w, x, y, z));
            case 3: return derivative(fun, wrt(z), at(w, x, y, z));
            default: throw std::invalid_argument("wrtidx must be 0, 1, 2 or 3");
        }
    }, py::arg("func"), py::arg("w"), py::arg("x"), py::arg("y"), py::arg("z"), py::arg("wrt"), py::return_value_policy::reference);
}

void export_derivative1(py::module& m) { exportDerivative1<std::function<real(real)>, real>(m, "derivative1"); }
void export_derivative2(py::module& m) { exportDerivative2<std::function<real(real, real)>, real>(m, "derivative2"); }
void export_derivative3(py::module& m) { exportDerivative3<std::function<real(real, real, real)>, real>(m, "derivative3"); }
void export_derivative4(py::module& m) { exportDerivative4<std::function<real(real, real, real, real)>, real>(m, "derivative4"); }