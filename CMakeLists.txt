cmake_minimum_required(VERSION 3.21)

# Suppress developer warnings from external packages
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS
    ON
    CACHE INTERNAL ""
)

# Set CMake policies to suppress warnings
if(POLICY CMP0148)
  cmake_policy(SET CMP0148 NEW)
endif()
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
if(POLICY CMP0177)
  cmake_policy(SET CMP0177 NEW)
endif()

project(LuPNT LANGUAGES CXX)

# Avoid conflicts with Homebrew and MacPorts on macOS
set(CMAKE_IGNORE_PREFIX_PATH "/opt/homebrew;/usr/local")

include(${CMAKE_CURRENT_LIST_DIR}/cmake/Tools.cmake)

# # needed to generate test target enable_testing()

message(STATUS "Configuring library...")
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/cpp ${CMAKE_BINARY_DIR}/lupnt)

message(STATUS "Configuring examples...")
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/cpp/examples ${CMAKE_BINARY_DIR}/examples)

message(STATUS "Configuring python bindings...")
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/python/bindings ${CMAKE_BINARY_DIR}/bindings)

message(STATUS "Configuring projects...")
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/projects ${CMAKE_BINARY_DIR}/projects)
# ${CMAKE_BINARY_DIR}/projects) add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../test
# ${CMAKE_BINARY_DIR}/test) add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../projects
# ${CMAKE_BINARY_DIR}/projects)
