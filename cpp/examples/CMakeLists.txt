cmake_minimum_required(VERSION 3.21)

project(lupnt_examples LANGUAGES CXX)

# ---- Dependencies ----
set(LUPNT_WITH_PYTHON
    OFF
    CACHE BOOL "Build with Python interpreter." FORCE
)

# Check if lupnt target already exists (e.g., from parent CMakeLists.txt)
include(${CMAKE_CURRENT_LIST_DIR}/../../cmake/CPM.cmake)
if(NOT TARGET lupnt::lupnt)
  CPMAddPackage(NAME lupnt SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)
endif()

# Add TBB ------ Make conda discoverable for find_package
list(PREPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
set(TBB_DIR "$ENV{CONDA_PREFIX}/lib/cmake/TBB")
find_package(TBB CONFIG REQUIRED COMPONENTS tbb tbbmalloc)

# ---- Create examples executable ----

file(
  GLOB_RECURSE examples
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  *.cc
)

foreach(example ${examples})
  get_filename_component(name ${example} NAME_WE)
  get_filename_component(dir ${example} DIRECTORY)

  add_executable(${name} ${example})
  # Link to the lupnt target (either existing or created by CPM)
  target_link_libraries(${name} PRIVATE lupnt::lupnt TBB::tbb TBB::tbbmalloc)
  target_link_directories(${name} PRIVATE "$ENV{CONDA_PREFIX}/lib")
  set_target_properties(${name} PROPERTIES CXX_STANDARD 20)
  set_target_properties(${name} PROPERTIES CXX_STANDARD_REQUIRED ON)
  set_target_properties(${name} PROPERTIES FOLDER "examples/${dir}")

  # if MACOS
  if(APPLE)
    set_target_properties(
      ${name}
      PROPERTIES
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH
        "@loader_path;\
@loader_path/../lib;\
@loader_path/../lupnt;\
@loader_path/../lupnt/lupnt/environment/plasma/fortran;\
$ENV{CONDA_PREFIX}/lib"
    )
  endif()

  # Fix macOS SDK compatibility issues with conda clang
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND APPLE)
    target_compile_definitions(${name} PRIVATE _LIBCPP_DISABLE_AVAILABILITY)
    target_compile_options(
      ${name} PRIVATE -Wno-error=availability -Wno-error=deprecated-non-prototype
    )
  endif()
endforeach()

# Create a custom target to build all examples
add_custom_target(all_examples)
foreach(example ${examples})
  get_filename_component(name ${example} NAME_WE)
  add_dependencies(all_examples ${name})
  set_target_properties(all_examples PROPERTIES FOLDER "examples")
endforeach()
