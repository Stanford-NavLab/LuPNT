cmake_minimum_required(VERSION 3.15)
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON)

project(LuPNTExamples LANGUAGES CXX)

# --- Import tools ----

include(${CMAKE_CURRENT_LIST_DIR}/../../cmake/tools.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/FindOpenMP.cmake)

# ---- Dependencies ----

include(${CMAKE_CURRENT_LIST_DIR}/../../cmake/CPM.cmake)
CPMAddPackage(NAME LuPNT SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/../..)

# ---- Create examples executable ----

file(
  GLOB_RECURSE examples
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  *.cc
)
# Remove example_ekf_elfo.cc from the list to compile
list(REMOVE_ITEM examples "example_ekf_elfo.cc")

find_package(OpenMP REQUIRED)

foreach(example ${examples})
  get_filename_component(name ${example} NAME_WE)
  get_filename_component(dir ${example} DIRECTORY)

  add_executable(${name} ${example})
  target_link_libraries(${name} PRIVATE LuPNT::LuPNT)
  set_target_properties(${name} PROPERTIES CXX_STANDARD 20)
  set_target_properties(${name} PROPERTIES FOLDER "examples/${dir}")
  set_target_properties(
    ${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/examples
  )
endforeach()

# Create a custom target to build all examples
add_custom_target(all_examples)
foreach(example ${examples})
  get_filename_component(name ${example} NAME_WE)
  add_dependencies(all_examples ${name})
  set_target_properties(all_examples PROPERTIES FOLDER "examples")
endforeach()
