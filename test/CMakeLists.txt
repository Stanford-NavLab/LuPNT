cmake_minimum_required(VERSION 3.15)
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON)

project(LuPNTTests LANGUAGES CXX)

# ---- Options ----

option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
option(TEST_INSTALLED_VERSION "Test the version found by find_package" OFF)

# --- Import tools ----

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/tools.cmake)

# ---- Dependencies ----
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/FindOpenMP.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/CPM.cmake)

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
  find_package(LuPNT REQUIRED)
else()
  CPMAddPackage(NAME LuPNT SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)
endif()

# ---- Create binary ----

set(ENV{CATCH_VERBOSITY} 2)
file(GLOB_RECURSE sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cpp/*.cc)
add_executable(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} Catch2::Catch2 LuPNT::LuPNT)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 20)

# ---- Add LuPNTTests ----

enable_testing()

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
message(STATUS "Catch2_SOURCE_DIR: ${Catch2_SOURCE_DIR}")
include(CTest)
include(Catch)
catch_discover_tests(${PROJECT_NAME})

# ---- code coverage ----

if(ENABLE_TEST_COVERAGE)
  target_compile_options(LuPNT PUBLIC -O0 -g -fprofile-arcs -ftest-coverage)
  target_link_options(LuPNT PUBLIC -fprofile-arcs -ftest-coverage)
endif()
