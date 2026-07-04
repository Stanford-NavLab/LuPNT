cmake_minimum_required(VERSION 3.21)

# Set CMake policies to suppress warnings
if(POLICY CMP0177)
  cmake_policy(SET CMP0177 NEW)
endif()

# Avoid conflicts with Homebrew and MacPorts on macOS
set(CMAKE_IGNORE_PREFIX_PATH "/opt/homebrew;/usr/local")

project(lupnt_python LANGUAGES CXX)

# ---- Dependencies ----
include(${CMAKE_CURRENT_LIST_DIR}/../../cmake/CPM.cmake)

# Add TBB ------ Make conda discoverable for find_package
list(PREPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
set(TBB_DIR "$ENV{CONDA_PREFIX}/lib/cmake/TBB")
find_package(TBB CONFIG REQUIRED COMPONENTS tbb tbbmalloc)

set(LUPNT_WITH_PYTHON
    OFF
    CACHE BOOL "Build with Python interpreter." FORCE
)

if(NOT TARGET lupnt::lupnt)
  CPMAddPackage(NAME lupnt SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/../../cpp)
endif()
# Use FindPython instead of deprecated FindPythonInterp/FindPythonLibs
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 REQUIRED)

# ---- Create standalone executable ----
file(
  GLOB_RECURSE cc_files
  RELATIVE ${CMAKE_CURRENT_LIST_DIR}
  *.cc
)

pybind11_add_module(_pylupnt MODULE ${cc_files})
target_link_libraries(_pylupnt PRIVATE lupnt::lupnt TBB::tbb TBB::tbbmalloc)

# Ensure the linker searches your env’s lib dir (fixes -ltbb)
target_link_directories(_pylupnt PRIVATE "$ENV{CONDA_PREFIX}/lib")

# if MACOS
if(APPLE)
  set_target_properties(
    _pylupnt PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE
                        INSTALL_RPATH "@loader_path;@loader_path/../lib;$ENV{CONDA_PREFIX}/lib"
  )
endif()
set_target_properties(_pylupnt PROPERTIES CXX_STANDARD 20)

# Disable sanitizers for Python module (Python can't import modules with sanitizers) Explicitly
# disable all sanitizers - this flag overrides global sanitizer settings
target_compile_options(_pylupnt PRIVATE -fno-sanitize=all)
if(NOT APPLE)
  target_link_options(_pylupnt PRIVATE -fno-sanitize=all)
endif()

# Fix macOS SDK compatibility issues with conda clang
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND APPLE)
  target_compile_definitions(_pylupnt PRIVATE _LIBCPP_DISABLE_AVAILABILITY)
  target_compile_options(
    _pylupnt PRIVATE -Wno-error=availability -Wno-error=deprecated-non-prototype
  )
endif()

# Fix linker issues
if(APPLE)
  # Set proper RPATH handling to avoid duplicates
  set_target_properties(
    _pylupnt
    PROPERTIES LINK_FLAGS "-undefined dynamic_lookup"
               INSTALL_RPATH "@loader_path;@loader_path/../lib"
               BUILD_WITH_INSTALL_RPATH TRUE
  )

  # Ensure consistent visibility settings
  target_compile_options(_pylupnt PRIVATE -fvisibility=hidden)
  target_compile_definitions(_pylupnt PRIVATE -DHAVE_HIDDEN_VISIBILITY)
endif(APPLE)

# Ensure consistent C++ standard library visibility
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(_pylupnt PRIVATE -fvisibility-inlines-hidden)
endif()

install(TARGETS _pylupnt DESTINATION pylupnt)

# ---- Utility for development build ----

add_custom_target(
  pylupnt-dev
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target _pylupnt
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:_pylupnt> ${CMAKE_CURRENT_LIST_DIR}/../pylupnt
  COMMAND "$ENV{CONDA_PREFIX}/bin/stubgen" -m _pylupnt -o ${CMAKE_CURRENT_LIST_DIR}/../pylupnt
  WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/../pylupnt
  COMMENT "Building _pylupnt and generating new stubs."
)

# Install the C++ shared library to a location discoverable by _pylupnt
install(TARGETS lupnt LIBRARY DESTINATION pylupnt/../lib)

# For development: copy the dylib after build
add_custom_command(
  TARGET _pylupnt
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_LIST_DIR}/../lib
  COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:lupnt>
          ${CMAKE_CURRENT_LIST_DIR}/../lib/
  COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:_pylupnt>
          ${CMAKE_CURRENT_LIST_DIR}/../pylupnt/
  COMMENT "Syncing _pylupnt extension module to python/pylupnt/"
)

if(TARGET iri2007_lib)
  add_custom_command(
    TARGET _pylupnt
    POST_BUILD
    # Also copy libiri2007_lib.dylib (SHARED Fortran lib) so Python finds it via @loader_path in
    # liblupnt.dylib's rpath.
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:iri2007_lib>
            ${CMAKE_CURRENT_LIST_DIR}/../lib/
  )
endif()
