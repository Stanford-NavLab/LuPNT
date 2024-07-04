find_package(Eigen3 REQUIRED)
find_package(Catch2 3 REQUIRED)

file(GLOB_RECURSE TEST_CPP_FILES "test_*.cc")
add_executable(lupnt-cpptests main.cpp ${TEST_CPP_FILES})
target_link_libraries(lupnt-cpptests lupnt::lupnt Catch2::Catch2WithMain)
set_target_properties(lupnt-cpptests PROPERTIES CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)

# Note: Disabled compilation of implicit conversion tests because they are
# stressing too much CI compilers in addition to the other tests. 
# Stanford Navigation and Autonomous Vehicles (NAV) Lab, 19.06.2023

# add_executable(lupnt-cpptests-implicit-conversion main.cpp ${TEST_CPP_FILES})
# target_compile_definitions(lupnt-cpptests-implicit-conversion PUBLIC LUPNT_ENABLE_IMPLICIT_CONVERSION=1)
# target_link_libraries(lupnt-cpptests-implicit-conversion lupnt::lupnt Catch2::Catch2WithMain)
# set_target_properties(lupnt-cpptests-implicit-conversion PROPERTIES CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)

# Add target tests that performs all C++ and Python tests
add_custom_target(tests
        COMMENT "Running C++ tests..."
        COMMAND $<TARGET_FILE:lupnt-cpptests>
        # COMMAND $<TARGET_FILE:lupnt-cpptests-implicit-conversion>
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

if(LUPNT_TEST_SANITIZE)
    include(CheckCXXCompilerFlag)
    set(CMAKE_REQUIRED_FLAGS "-fsanitize=address")  # Also needs to be a link flag for test to work
    check_cxx_compiler_flag(-fsanitize=address HAVE_ASAN)
    unset(CMAKE_REQUIRED_FLAGS)
    if(HAVE_ASAN)
        target_compile_options(lupnt-cpptests PRIVATE "-fsanitize=address")
        # target_compile_options(lupnt-cpptests-implicit-conversion PRIVATE "-fsanitize=address")
        target_link_options(lupnt-cpptests PRIVATE "-fsanitize=address")
        # target_link_options(lupnt-cpptests-implicit-conversion PRIVATE "-fsanitize=address")
    endif()
endif()
