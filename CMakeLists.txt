# Project setup
cmake_minimum_required(VERSION 3.15..3.22)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 20) # C++20...
set(CMAKE_CXX_STANDARD_REQUIRED ON) # ...is required...
set(CMAKE_CXX_EXTENSIONS OFF) # ...without compiler extensions like gnu++11
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
add_definitions(-D_USE_MATH_DEFINES)

project(LuPNT LANGUAGES CXX)

# set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")

# guard against in-source buildsnvi 
if(${LuPNT_SOURCE_DIR} STREQUAL ${LuPNT_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

# build directories
add_subdirectory(third_party)
add_subdirectory(lupnt)
add_subdirectory(examples/cpp) 
add_subdirectory(bindings)

option (BUILD_TESTING "Build the testing tree." ON)
# # Only build tests if we are the top-level project
# # Allows this to be used by super projects with `add_subdirectory`
if (BUILD_TESTING AND (PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
  enable_testing()
  add_subdirectory(test/cpp)
endif()


# Codecov
set(CODECOV_TOKEN $ENV{CODECOV_TOKEN})
set(ENV{PATH} "$ENV{PATH}:/usr/local/opt/llvm/bin")

add_custom_target(codecov
  # Run tests
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target tests
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test/cpp

  # Generate coverage data
  COMMAND /usr/local/opt/llvm/bin/llvm-profdata merge -sparse default.profraw -o default.profdata
  COMMAND /usr/local/opt/llvm/bin/llvm-cov show ${CMAKE_BINARY_DIR}/lupnt/liblupnt.a -instr-profile=default.profdata > coverage.txt

  # Download and execute Codecov uploader
  COMMAND curl -Os https://uploader.codecov.io/latest/macos/codecov
  COMMAND chmod +x codecov
  COMMAND ./codecov -t ${CODECOV_TOKEN} -f coverage.txt

  COMMENT "Running tests, generating coverage report and uploading to Codecov"
)
