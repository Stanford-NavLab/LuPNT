cmake_minimum_required(VERSION 3.14...3.22)
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON)

project(LuPNTStandalone LANGUAGES CXX)

# --- Import tools ----

include(${CMAKE_CURRENT_LIST_DIR}/../cmake/tools.cmake)

# ---- Dependencies ----

include(${CMAKE_CURRENT_LIST_DIR}/../cmake/CPM.cmake)

CPMAddPackage(
  GITHUB_REPOSITORY jarro2783/cxxopts
  VERSION 3.0.0
  OPTIONS "CXXOPTS_BUILD_EXAMPLES NO" "CXXOPTS_BUILD_TESTS NO" "CXXOPTS_ENABLE_INSTALL YES"
)

CPMAddPackage(NAME LuPNT SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)

# ---- Create standalone executable ----

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)

add_executable(${PROJECT_NAME} ${sources})

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 20 OUTPUT_NAME "LuPNT")

target_link_libraries(${PROJECT_NAME} LuPNT::LuPNT cxxopts)
