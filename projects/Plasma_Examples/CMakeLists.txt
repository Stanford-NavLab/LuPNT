cmake_minimum_required(VERSION 3.15)
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON)

get_filename_component(project_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${project_name} LANGUAGES CXX)

include(${CMAKE_CURRENT_LIST_DIR}/../../cmake/CPM.cmake)

list(PREPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
set(TBB_DIR "$ENV{CONDA_PREFIX}/lib/cmake/TBB")
find_package(TBB CONFIG REQUIRED COMPONENTS tbb tbbmalloc)
