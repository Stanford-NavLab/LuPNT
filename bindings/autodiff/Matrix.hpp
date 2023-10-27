/**
 * @file matrix.hpp
 * @author Keidai Iiyama 
 * @brief Utility functions for converting Eigen::MatrixXd to numpy array
 * @version 0.1
 * @date 2023-07-12
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#pragma once

#include "../pybind11.hxx"
#include <Eigen/Core>
#include <iostream>

/**
 * @brief Convert Eigen::MatrixXd to numpy array
 *        Reference: https://stackoverflow.com/questions/44659924/returning-numpy-arrays-via-pybind11
 * @param J 
 * @return py::array_t<double> 
 */
py::array_t<double> EigenMatrixToPythonArray(const Eigen::MatrixXd& J){
    // directly return a numpy array
    double *foo = new double[J.size()];
    for (size_t i = 0; i < J.size(); i++) {
        foo[i] = (double) *(J.data() + i);
    }

    // Create a Python object that will free the allocated
    // memory when destroyed:
    py::capsule free_when_done(foo, [](void *f) {
        double *foo = reinterpret_cast<double *>(f);
        // std::cerr << "Element [0] = " << foo[0] << "\n";
        // std::cerr << "freeing memory @ " << f << "\n";
        delete[] foo;
    });

    auto result = py::array_t<double>(
        {J.rows(), J.cols()},   /* Pointer to data (nullptr -> ask NumPy to allocate!) */
        {sizeof(double), sizeof(double)*J.rows()}, /* Strides (in bytes) for each index */
        foo,
        free_when_done /* Pointer to parent */
    );

    return result;
}