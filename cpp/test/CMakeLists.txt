cmake_minimum_required(VERSION 3.21)

project(lupnt_tests LANGUAGES CXX)

# ---- Options ----

option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
option(TEST_INSTALLED_VERSION "Test the version found by find_package" OFF)

foreach(_prefix "$ENV{HOME}/anaconda3" "$ENV{LUPNT_PYTHON_PATH}")
  if(IS_DIRECTORY "${_prefix}")
    list(APPEND CMAKE_IGNORE_PREFIX_PATH "${_prefix}")
    list(APPEND CMAKE_SYSTEM_IGNORE_PREFIX_PATH "${_prefix}")
  endif()
endforeach()

# --- Import tools ----
# NOTE: this file lives at `cpp/test/`, two directories below the project root, so the shared CMake
# helpers (which live in `<root>/cmake/`, see `cpp/CMakeLists.txt`) are reached via `../../cmake/`,
# not `../cmake/` (which would resolve to the nonexistent `cpp/cmake/`).

set(LUPNT_TOOLS_CMAKE "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/Tools.cmake")
if(EXISTS "${LUPNT_TOOLS_CMAKE}")
  include("${LUPNT_TOOLS_CMAKE}")
endif()

# ---- Dependencies ----
# `FindOpenMP.cmake`/`Python.cmake` don't exist anywhere in the repo (and aren't needed here):
# OpenMP is located via the standard CMake `FindOpenMP` module by `find_package(OpenMP REQUIRED
# COMPONENTS CXX)` in `cpp/CMakeLists.txt`, and is pulled in transitively when `lupnt::lupnt` is
# built below via `CPMAddPackage`.
include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/CPM.cmake)

CPMAddPackage("gh:TheLartians/Format.cmake@1.7.3")
set_target_properties(check-cmake-format PROPERTIES FOLDER format)
set_target_properties(check-clang-format PROPERTIES FOLDER format)
set_target_properties(check-format PROPERTIES FOLDER format)
set_target_properties(fix-cmake-format PROPERTIES FOLDER format)
set_target_properties(fix-clang-format PROPERTIES FOLDER format)
set_target_properties(fix-format PROPERTIES FOLDER format)
set_target_properties(clang-format PROPERTIES FOLDER format)
set_target_properties(cmake-format PROPERTIES FOLDER format)
set_target_properties(format PROPERTIES FOLDER format)

CPMAddPackage(
  NAME Catch2
  GITHUB_REPOSITORY catchorg/Catch2
  VERSION 3.6.0
  # OPTIONS "CMAKE_BUILD_TYPE=Release"
)
set_target_properties(Catch2 PROPERTIES FOLDER thirdparty)
set_target_properties(Catch2WithMain PROPERTIES FOLDER thirdparty)

if(TEST_INSTALLED_VERSION)
  find_package(lupnt REQUIRED)
else()
  CPMAddPackage(NAME lupnt SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)
endif()

# ---- Create binary ----

set(ENV{CATCH_VERBOSITY} 2)
# NOTE: the test sources live directly under `cpp/test/` (and its subdirectories, e.g. `agents/`,
# `interfaces/`, ...), not under a nested `cpp/test/cpp/` -- the old `.../cpp/*.cc` pattern matched
# nothing.
file(GLOB_RECURSE sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)
add_executable(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} Catch2::Catch2 lupnt::lupnt)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 20)
if(TARGET iri2007_lib)
  set_property(
    TARGET ${PROJECT_NAME}
    APPEND
    PROPERTY BUILD_RPATH "$<TARGET_FILE_DIR:iri2007_lib>"
  )
endif()

# Fix macOS SDK compatibility issues with conda clang
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND APPLE)
  target_compile_definitions(${PROJECT_NAME} PRIVATE _LIBCPP_DISABLE_AVAILABILITY)
  target_compile_options(
    ${PROJECT_NAME} PRIVATE -Wno-error=availability -Wno-error=deprecated-non-prototype
  )
endif()

# ---- Add lupnt_tests ----

enable_testing()

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
message(STATUS "Catch2_SOURCE_DIR: ${Catch2_SOURCE_DIR}")
include(CTest)
include(Catch)
catch_discover_tests(${PROJECT_NAME})

# ---- code coverage ----

if(ENABLE_TEST_COVERAGE)
  target_compile_options(lupnt PUBLIC -O0 -g -fprofile-arcs -ftest-coverage)
  target_link_options(lupnt PUBLIC -fprofile-arcs -ftest-coverage)
endif()
