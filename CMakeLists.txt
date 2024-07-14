# The minimum required cmake version
cmake_minimum_required(VERSION 3.20)

# Add cmake modules of this project to the module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Use ccache to speed up repeated compilations
include(CCache)

# Name and details of the project
project(lupnt VERSION 1.1.2 LANGUAGES CXX)

# Include the cmake variables with values for installation directories
include(GNUInstallDirs)

# Ensure proper configuration if in a conda environment
include(CondaAware)

# Define build options
option(LUPNT_BUILD_TESTS "Enable the compilation of the test files." ON)
option(LUPNT_BUILD_PYTHON "Enable the compilation of the python bindings." ON)
option(LUPNT_BUILD_EXAMPLES "Enable the compilation of the example files." ON)
option(LUPNT_BUILD_DOCS "Enable the build of the documentation and website." ON)

# lupnt requires a c++23 enabled compiler
set(CMAKE_CXX_STANDARD 23)          # ensure cmake instructs compiler to use C++23
set(CMAKE_CXX_STANDARD_REQUIRED ON) # ensure the C++ standard given before is actually used
set(CMAKE_CXX_EXTENSIONS OFF)       # avoid compile flags of the type -std=gnu++1z added by cmake

# Set compiler warning flags
# if(MSVC)
#     add_compile_options(/W4 /WX) # Use /W4 for level 4 warnings, /WX to treat warnings as errors
# else()
#     add_compile_options(-Wall -Wextra -Wpedantic -Werror) # For GCC and Clang
# endif()

# OpenMP
set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
set(OpenMP_CXX_LIB_NAMES "omp")

# Build the library
add_subdirectory(lupnt)

# Build the tests
if(LUPNT_BUILD_TESTS)
    add_subdirectory(tests/cpp)
endif()

# Build the python wrappers
if(LUPNT_BUILD_PYTHON)
    find_package(Python COMPONENTS Interpreter Development)
    find_package(pybind11 REQUIRED)
    add_subdirectory(python)
endif()

# # # Build the examples
if(LUPNT_BUILD_EXAMPLES)
    add_subdirectory(examples/cpp)
endif()

# # Build the docs
# if(LUPNT_BUILD_DOCS)
#     add_subdirectory(docs)
# endif()

# Install the cmake config files that permit users to use find_package(lupnt)
include(lupntInstallCMakeConfigFiles)
