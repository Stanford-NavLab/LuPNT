################################
# Testing
################################
if (APPLE)
  add_definitions(-DGTEST_USE_OWN_TR1_TUPLE)
  add_definitions(-D__GLIBCXX__)
endif (APPLE)

#   # This adds another subdirectory, which has 'project(gtest)'.
#   add_subdirectory(lib/gtest-1.6.0)

#   enable_testing()

#   # Include the gtest library. gtest_SOURCE_DIR is available due to
#   # 'project(gtest)' above.
#   include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

##############
# Unit Tests
##############
enable_testing()

file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cc)
set(SOURCES ${TEST_SOURCES})
add_executable(runUnitTests  ${TEST_SOURCES})

target_include_directories(runUnitTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../lupnt)

# Standard linking to gtest stuff.
target_link_libraries(runUnitTests GTest::gtest_main)

# Extra linking for the project.
target_link_libraries(runUnitTests lupnt)

# This is so you can do 'make test' to see all your tests run, instead of
# manually running the executable runUnitTests to see those specific tests.
# find_package(GTest REQUIRED)
include(GoogleTest)
gtest_discover_tests(runUnitTests)

# You can also omit NAME and COMMAND. The second argument could be some other
# test executable.
#add_test(that-other-test-I-made runUnitTests)