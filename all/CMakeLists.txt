# this script adds all subprojects to a single build to allow IDEs understand the full project
# structure.

cmake_minimum_required(VERSION 3.15)
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON)

project(BuildAll LANGUAGES CXX)

include(${CMAKE_CURRENT_LIST_DIR}/../cmake/tools.cmake)

# needed to generate test target
enable_testing()

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../test ${CMAKE_BINARY_DIR}/test)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../source/python ${CMAKE_BINARY_DIR}/python)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../examples/cpp ${CMAKE_BINARY_DIR}/examples)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../projects ${CMAKE_BINARY_DIR}/projects)
