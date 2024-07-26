cmake_minimum_required(VERSION 3.14...3.22)
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON)

# ---- Project ----

project(
  LuPNT
  VERSION 1.0
  LANGUAGES CXX
)

# ---- Include guards ----

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

# ---- Add dependencies ----

include(cmake/FindOpenMP.cmake)
include(cmake/CPM.cmake)
add_subdirectory(thirdparty)

# ---- Add source files ----

# Note: globbing sources is considered bad practice as CMake's generators may not detect new files
# automatically. Keep that in mind when changing files, or explicitly mention them here.
file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/source/*.cc")

# ---- Create library ----

# Note: for header-only libraries change all PUBLIC flags to INTERFACE and create an interface
# target: add_library(${PROJECT_NAME} INTERFACE)
add_library(${PROJECT_NAME} ${headers} ${sources})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)

# being a cross-platform target, we enforce standards conformance on MSVC
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/permissive->")

function(target_link_libraries_system target)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM PUBLIC ${lib_include_dirs})
    target_link_libraries(${target} PUBLIC SYSTEM ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

# Link dependencies
target_link_libraries_system(
  ${PROJECT_NAME} fmt::fmt Eigen3::Eigen matplot autodiff::autodiff cspice
)

target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
target_include_directories(${PROJECT_NAME} PUBLIC ${OpenMP_INCLUDE_DIRS})

target_include_directories(
  ${PROJECT_NAME}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
         $<BUILD_INTERFACE:${cspice_SOURCE_DIR}/include>
)

# enable compiler warnings if(NOT TEST_INSTALLED_VERSION) if(CMAKE_CXX_COMPILER_ID MATCHES "Clang"
# OR CMAKE_CXX_COMPILER_ID MATCHES "GNU") target_compile_options(${PROJECT_NAME} PRIVATE -Wall
# -Wpedantic -Wextra -Werror) elseif(MSVC) # target_compile_options(${PROJECT_NAME} PRIVATE /W4 /WX)
# target_compile_definitions(${PROJECT_NAME} PRIVATE DOCTEST_CONFIG_USE_STD_HEADERS) endif() endif()

# ---- Create an installable target ----
# this allows users to install and find the library via `find_package()`.

# the location where the project's version header will be placed should match the project's regular
# header paths
string(TOLOWER ${PROJECT_NAME}/version.h VERSION_HEADER_LOCATION)

packageProject(
  # the name of the target to export
  NAME ${PROJECT_NAME}
  # the version of the target to export
  VERSION ${PROJECT_VERSION}
  # a temporary directory to create the config files
  BINARY_DIR ${PROJECT_BINARY_DIR}
  # location of the target's public headers
  INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include
  # should match the target's INSTALL_INTERFACE include directory
  INCLUDE_DESTINATION
    include/${PROJECT_NAME}-${PROJECT_VERSION}
    # (optional) option to install only header files with matching pattern
    INCLUDE_HEADER_PATTERN "*.h"
    # semicolon separated list of the project's dependencies
  DEPENDENCIES ${LUPNT_DEPENDENCIES}
  # (optional) create a header containing the version info Note: that the path to headers should be
  # lowercase
  VERSION_HEADER "${VERSION_HEADER_LOCATION}"
  # (optional) install your library with a namespace (Note: do NOT add extra '::')
  NAMESPACE ${PROJECT_NAME}
  # (optional) define the project's version compatibility, defaults to `AnyNewerVersion` supported
  # values: `AnyNewerVersion|SameMajorVersion|SameMinorVersion|ExactVersion`
  COMPATIBILITY
    SameMajorVersion
    # (optional) option to disable the versioning of install destinations
    DISABLE_VERSION_SUFFIX
    YES
    # (optional) option to ignore target architecture for package resolution defaults to YES for
    # header only (i.e. INTERFACE) libraries
    ARCH_INDEPENDENT
    YES
    # (optional) option to generate CPack variables
    CPACK
    YES
)
