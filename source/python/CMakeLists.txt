cmake_minimum_required(VERSION 3.15)
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON)

project(LuPNTPython LANGUAGES CXX)

# --- Import tools ----

include(${CMAKE_CURRENT_LIST_DIR}/../../cmake/tools.cmake)

# ---- Dependencies ----
include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/FindOpenMP.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/../../cmake/CPM.cmake)
CPMAddPackage("gh:pybind/pybind11@2.13.1")
CPMAddPackage(NAME LuPNT SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/../..)

# ---- Create standalone executable ----

# file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/bindings/*.cc)
pybind11_add_module(_pylupnt bindings/py_main.cc)
target_link_libraries(_pylupnt PUBLIC LuPNT::LuPNT)
set_target_properties(_pylupnt PROPERTIES CXX_STANDARD 20)
if(APPLE)
  set_target_properties(_pylupnt PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(APPLE)
install(TARGETS _pylupnt DESTINATION pylupnt)

# ---- Utility for development build ----

add_custom_target(
  pylupnt-dev
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target _pylupnt
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:_pylupnt> .
  COMMAND stubgen -m _pylupnt -o .
  WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/pylupnt
  COMMENT
    "Building _pylupnt, removing old stubs, generating new stubs, and copying _pylupnt to current directory"
)
