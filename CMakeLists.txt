cmake_minimum_required(VERSION 3.15..3.22)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 20) # C++20...
set(CMAKE_CXX_STANDARD_REQUIRED ON) # ...is required...
set(CMAKE_CXX_EXTENSIONS OFF) # ...without compiler extensions like gnu++11
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
add_definitions(-D_USE_MATH_DEFINES)

# for conda
# https://conda-forge.org/docs/maintainer/knowledge_base.html#newer-c-features-with-old-sdk
add_definitions(-D_LIBCPP_DISABLE_AVAILABILITY)

# # Set our project name
project(LPT
  LANGUAGES CXX)

# guard against in-source builds
if(${LPT_SOURCE_DIR} STREQUAL ${LPT_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# GoogleTest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
# add_library(GTest::GTest INTERFACE IMPORTED)
# target_link_libraries(GTest::GTest INTERFACE gtest_main)

# Libraries:
include(cmake/HandleBoost.cmake) # Boost
include(cmake/HandleEigen.cmake) # Eigen3

# build directories
add_subdirectory(lupnt) # src directory
if (BUILD_EXAMPLES)
  add_subdirectory(examples/cpp) # src directory
endif()
add_subdirectory(test/cpp) # test directory


# PYBIND11 ---- check
add_subdirectory(pybind11)
# FetchContent_Declare(
#   pybind11
#   GIT_REPOSITORY https://github.com/pybind/pybind11.git
#   GIT_TAG        v2.6.2
# )
# FetchContent_MakeAvailable(pybind11)
add_subdirectory(bindings)

# install python libraries
# execute_process(
#   COMMAND "pip3 install ."
# )
