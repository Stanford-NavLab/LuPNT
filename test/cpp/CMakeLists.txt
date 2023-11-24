# GoogleTest
# include(FetchContent)
# FetchContent_Declare(
#   googletest
#   URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
#   DOWNLOAD_EXTRACT_TIMESTAMP TRUE
# )

# # For Windows: Prevent overriding the parent project's compiler/linker settings
# set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# FetchContent_MakeAvailable(googletest)
# include(GoogleTest)

# # add_library(GTest::GTest INTERFACE IMPORTED)
# # target_link_libraries(GTest::GTest INTERFACE gtest_main)

# if (APPLE)
#   add_definitions(-DGTEST_USE_OWN_TR1_TUPLE)
#   add_definitions(-D__GLIBCXX__)
# endif (APPLE)

# file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cc)
# # add_executable(tests ${TEST_SOURCES})

# target_include_directories(tests PRIVATE ${CMAKE_SOURCE_DIR}/lupnt)

# # Standard linking to gtest stuff.
# target_link_libraries(tests GTest::gtest_main)

# # Extra linking for the project.
# target_link_libraries(tests lupnt)

# gtest_discover_tests(tests
#   PROPERTIES ENVIRONMENT "LUPNT_DATA_PATH=${LUPNT_DATA_PATH}"
# )

# message(STATUS "LUPNT_DATA_PATH: ${LUPNT_DATA_PATH}")

Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.4.0 # or a later release
)

FetchContent_MakeAvailable(Catch2)

file(GLOB_RECURSE test_sources LIST_DIRECTORIES false *.cc)
add_executable(tests ${test_sources})
target_include_directories(tests PRIVATE ${CMAKE_SOURCE_DIR}/lupnt)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain lupnt)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
catch_discover_tests(tests
  PROPERTIES ENVIRONMENT "LUPNT_DATA_PATH=${LUPNT_DATA_PATH}"
)
message(STATUS "LUPNT_DATA_PATH: ${LUPNT_DATA_PATH}")