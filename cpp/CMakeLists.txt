cmake_minimum_required(VERSION 3.21)

# ---- Project ----

project(lupnt LANGUAGES CXX Fortran)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Enable ccache if available for faster rebuilds
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
  message(STATUS "ccache found: ${CCACHE_PROGRAM}")
endif()

option(LUPNT_WITH_PYTHON "Build with Python interpreter." OFF)

foreach(_prefix "$ENV{HOME}/anaconda3" "$ENV{LUPNT_PYTHON_PATH}")
  if(IS_DIRECTORY "${_prefix}")
    list(APPEND CMAKE_IGNORE_PREFIX_PATH "${_prefix}")
    list(APPEND CMAKE_SYSTEM_IGNORE_PREFIX_PATH "${_prefix}")
  endif()
endforeach()

# ---- Include guards ----

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

# ---- Add dependencies ----
include(${CMAKE_CURRENT_LIST_DIR}/../cmake/CondaAware.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/../cmake/CPM.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/../cmake/WriteHeader.cmake)

# ---- Fortran compiler flags for plasma/pecsim module ----
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-form -std=legacy -fno-automatic")
endif()
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/lupnt/environment/plasma/fortran)

# ---- Add source files ----

# Note: globbing sources is considered bad practice as CMake's generators may not detect new files
# automatically. Keep that in mind when changing files, or explicitly mention them here.
file(GLOB_RECURSE headers CONFIGURE_DEPENDS ${CMAKE_CURRENT_LIST_DIR}/lupnt/*.h)
file(GLOB_RECURSE sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_LIST_DIR}/lupnt/*.cc)
write_header(${CMAKE_CURRENT_LIST_DIR}/lupnt)

# message(STATUS "Headers:") foreach(h ${headers}) file(RELATIVE_PATH rel_h
# "${CMAKE_CURRENT_LIST_DIR}/lupnt" "${h}") message(STATUS "  ${rel_h}") endforeach()

# message(STATUS "Sources:") foreach(s ${sources}) file(RELATIVE_PATH rel_s
# "${CMAKE_CURRENT_LIST_DIR}/lupnt" "${s}") message(STATUS "  ${rel_s}") endforeach()

# ---- Set project version ----

file(READ "${CMAKE_CURRENT_LIST_DIR}/lupnt/version.txt" version_contents)
string(REGEX MATCH "MAJOR +([0-9]+)" _ ${version_contents})
set(PROJECT_VERSION_MAJOR ${CMAKE_MATCH_1})
string(REGEX MATCH "MINOR +([0-9]+)" _ ${version_contents})
set(PROJECT_VERSION_MINOR ${CMAKE_MATCH_1})
string(REGEX MATCH "PATCH +([0-9]+)" _ ${version_contents})
set(PROJECT_VERSION_PATCH ${CMAKE_MATCH_1})
set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")

# ---- Create library ----

# Note: for header-only libraries change all PUBLIC flags to INTERFACE and create an interface
# target: add_library(${PROJECT_NAME} INTERFACE)
add_library(${PROJECT_NAME} SHARED ${headers} ${sources})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 20)

# Enable Position Independent Code for linking into shared libraries
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Embed @loader_path in RPATH so that libiri2007_lib.dylib (SHARED) is found next to liblupnt.dylib
# when deployed to python/lib/
if(APPLE)
  set_target_properties(
    ${PROJECT_NAME} PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE INSTALL_RPATH
                                                             "@loader_path;$ENV{CONDA_PREFIX}/lib"
  )
endif()

# being a cross-platform target, we enforce standards conformance on MSVC
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/permissive->")

# Set compile definitions BEFORE precompiling headers
target_compile_definitions(
  ${PROJECT_NAME}
  PUBLIC MAGIC_ENUM_RANGE_MIN=0 MAGIC_ENUM_RANGE_MAX=1024 AUTODIFF_ENABLE_IMPLICIT_CONVERSION
         AUTODIFF_ENABLE_IMPLICIT_CONVERSION_REAL OPENCV_EIGEN_TENSOR_SUPPORT
)

# Enable precompiled headers for faster compilation
target_precompile_headers(
  ${PROJECT_NAME}
  PUBLIC
  <Eigen/Core>
  <Eigen/Dense>
  <autodiff/forward/real.hpp>
  <yaml-cpp/yaml.h>
  <nlohmann/json.hpp>
  <magic_enum/magic_enum.hpp>
)

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/thirdparty)

# Link Fortran plasma libraries
target_link_libraries(${PROJECT_NAME} PRIVATE gcpm_lib iri2007_lib iri2020_lib xform_lib igrf_lib)

# iri2007_lib is SHARED (to avoid duplicate Fortran symbols with iri2020_lib). Add it to the lupnt
# install/export so CMake's export set validation passes.
install(
  TARGETS iri2007_lib
  EXPORT lupntTargets
  LIBRARY DESTINATION lib
)

function(target_link_libraries_system)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${lib_include_dirs})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

# Link dependencies
target_link_libraries_system(matplot autodiff magic_enum::magic_enum cspice)

find_package(Eigen3 REQUIRED)
find_package(fmt REQUIRED)
# Use FindPython instead of deprecated FindPythonInterp/FindPythonLibs
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(spdlog REQUIRED)
find_package(OpenCV REQUIRED)
find_package(GDAL REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)

target_link_libraries(
  ${PROJECT_NAME}
  PUBLIC Eigen3::Eigen
         fmt::fmt
         pybind11::pybind11
         yaml-cpp::yaml-cpp
         spdlog::spdlog
         ${OpenCV_LIBS}
         GDAL::GDAL
         OpenMP::OpenMP_CXX
)

target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
                         $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)

message(STATUS "LUPNT_WITH_PYTHON: ${LUPNT_WITH_PYTHON}")
if(LUPNT_WITH_PYTHON)
  find_package(
    Python3
    COMPONENTS Interpreter Development
    REQUIRED
  )
  target_compile_definitions(${PROJECT_NAME} PRIVATE LUPNT_WITH_PYTHON)
endif()

# enable compiler warnings
function(enable_warnings target)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wno-sign-compare)
  elseif(MSVC)
    target_compile_options(${target} PRIVATE /W4)
    target_compile_definitions(${target} PRIVATE DOCTEST_CONFIG_USE_STD_HEADERS)
  endif()
endfunction()

enable_warnings(${PROJECT_NAME})

# Fix macOS SDK compatibility issues with conda clang
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND APPLE)
  target_compile_definitions(${PROJECT_NAME} PRIVATE _LIBCPP_DISABLE_AVAILABILITY)
  target_compile_options(
    ${PROJECT_NAME} PRIVATE -Wno-error=availability -Wno-error=deprecated-non-prototype
  )
endif()

# # Fix visibility settings to prevent weak symbol warnings if(CMAKE_CXX_COMPILER_ID MATCHES
# "Clang") target_compile_options(${PROJECT_NAME} PRIVATE -fvisibility=hidden
# -fvisibility-inlines-hidden) target_compile_definitions(${PROJECT_NAME} PRIVATE
# -DHAVE_HIDDEN_VISIBILITY) endif()

# ---- Create an installable target ----
# this allows users to install and find the library via `find_package()`.

# the location where the project's version header will be placed should match the project's regular
# header paths
string(TOLOWER ${PROJECT_NAME}/version.h VERSION_HEADER_LOCATION)

packageProject(
  # the name of the target to export
  NAME ${PROJECT_NAME}
  # the version of the target to export
  VERSION ${PROJECT_VERSION}
  # a temporary directory to create the config files
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}
  # location of the target's public headers
  INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/lupnt
  # should match the target's INSTALL_INTERFACE include directory
  INCLUDE_DESTINATION
    include/${PROJECT_NAME}-${PROJECT_VERSION}
    # (optional) option to install only header files with matching pattern
    INCLUDE_HEADER_PATTERN "*.h"
    # semicolon separated list of the project's dependencies
  DEPENDENCIES ${LUPNT_DEPENDENCIES}
  # (optional) create a header containing the version info Note: that the path to headers should be
  # lowercase
  VERSION_HEADER "${VERSION_HEADER_LOCATION}"
  # (optional) install your library with a namespace (Note: do NOT add extra '::')
  NAMESPACE ${PROJECT_NAME}
  # (optional) define the project's version compatibility, defaults to `AnyNewerVersion` supported
  # values: `AnyNewerVersion|SameMajorVersion|SameMinorVersion|ExactVersion`
  COMPATIBILITY
    SameMajorVersion
    # (optional) option to disable the versioning of install destinations
    DISABLE_VERSION_SUFFIX
    YES
    # (optional) option to ignore target architecture for package resolution defaults to YES for
    # header only (i.e. INTERFACE) libraries
    ARCH_INDEPENDENT
    YES
    # (optional) option to generate CPack variables
    CPACK
    YES
)
