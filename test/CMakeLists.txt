cmake_minimum_required(VERSION 3.14...3.22)
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON)

project(LuPNTTests LANGUAGES CXX)

# ---- Options ----

option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)
option(TEST_INSTALLED_VERSION "Test the version found by find_package" OFF)

# --- Import tools ----

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/tools.cmake)

# ---- Dependencies ----

include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/CPM.cmake)

CPMAddPackage("gh:doctest/doctest@2.4.9")
CPMAddPackage("gh:TheLartians/Format.cmake@1.7.3")

if(TEST_INSTALLED_VERSION)
  find_package(LuPNT REQUIRED)
else()
  CPMAddPackage(NAME LuPNT SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)
endif()

# ---- Create binary ----

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
add_executable(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} doctest::doctest LuPNT::LuPNT)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 20)

# ---- Add LuPNTTests ----

enable_testing()

# Note: doctest and similar testing frameworks can automatically configure CMake tests. For other
# testing frameworks add the tests target instead: add_test(NAME ${PROJECT_NAME} COMMAND
# ${PROJECT_NAME})

include(${doctest_SOURCE_DIR}/scripts/cmake/doctest.cmake)
doctest_discover_tests(${PROJECT_NAME})

# ---- code coverage ----

if(ENABLE_TEST_COVERAGE)
  target_compile_options(LuPNT PUBLIC -O0 -g -fprofile-arcs -ftest-coverage)
  target_link_options(LuPNT PUBLIC -fprofile-arcs -ftest-coverage)
endif()
