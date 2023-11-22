project(lupnt
  LANGUAGES CXX)

# We split the library in to separate subfolders, each containing
# tests, timing, and an optional convenience library.
# The following variable is the master list of subdirs to add
set(LPT_subdirs
  core
  dynamics
  numerics
  agents
  physics
  measurements
)

set(LPT_srcs)

# Build 3rdparty separately
message(STATUS "Building 3rdparty")
add_subdirectory(3rdparty)

# Add the "lupnt" directory to the include search paths
include_directories("${CMAKE_SOURCE_DIR}/lupnt")

set(3rdparty_srcs
  ${eigen_headers} # Set by 3rdparty/CMakeLists.txt
  ${autodiff_headers}
  ${cspice_headers}
  ${libInterpolate_headers}
)

# Add the parent directory to the header search path
# include_directories(${CMAKE_CURRENT_SOURCE_DIR})
# foreach(subdir ${LPT_subdirs})
#   include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${subdir})
# endforeach(subdir)

# assemble core libaries
foreach(subdir ${LPT_subdirs})
  # Build convenience libraries
  file(GLOB_RECURSE subdir_srcs "${subdir}/*.cc" "${subdir}/*.h") # Include header files so they show up in Visual Studio
  set(${subdir}_srcs ${subdir_srcs})

  # Build local library and tests
  message(STATUS "Building ${subdir}")
  add_subdirectory(${subdir})
endforeach(subdir)

# To add additional sources to gtsam when building the full library (static or shared)
# append the subfolder with _srcs appended to the end to this list
set(LPT_srcs ${3rdparty_srcs})

foreach(subdir ${LPT_subdirs})
  list(APPEND LPT_srcs ${${subdir}_srcs})
endforeach(subdir)

# build shared and static versions of the library
message(STATUS "Building lupnt - shared: ${BUILD_SHARED_LIBS}")

# BUILD_SHARED_LIBS automatically defines static/shared libs:
add_library(lupnt ${LPT_srcs})

# add cspice
target_link_libraries(lupnt PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/cspice/cspice.a")

# add Boost
target_link_libraries(lupnt PUBLIC ${Boost_LIBRARIES})

# include libraries
target_include_directories(lupnt BEFORE PUBLIC # main gtsam includes:
  $<BUILD_INTERFACE:${LPT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/>
)

# 3rdparty libraries: use the "system" flag so they are included via "-isystem"
# and warnings (and warnings-considered-errors) in those headers are not
# reported as warnings/errors in our targets:
target_include_directories(lupnt SYSTEM BEFORE PUBLIC

  # autodiff
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/autodiff>
  $<INSTALL_INTERFACE:include/lupnt/3rdparty/autodiff>

  # Eigen
  $<BUILD_INTERFACE:${LPT_EIGEN_INCLUDE_FOR_BUILD}>
  $<INSTALL_INTERFACE:${LPT_EIGEN_INCLUDE_FOR_INSTALL}>

  # cspice
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/cspice>
  $<INSTALL_INTERFACE:include/lupnt/3rdparty/cspice>

  # libInterpolate
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libInterpolate>
  $<INSTALL_INTERFACE:include/lupnt/3rdparty/libInterpolate>
)

if(WIN32) # Add 'lib' prefix to static library to avoid filename collision with shared library
  if(NOT BUILD_SHARED_LIBS)
    set_target_properties(lupnt PROPERTIES
      PREFIX "lib"
      COMPILE_DEFINITIONS LPT_IMPORT_STATIC)
  else()
    set_target_properties(lupnt PROPERTIES
      PREFIX ""
      DEFINE_SYMBOL LPT_EXPORTS
      RUNTIME_OUTPUT_DIRECTORY "${LPT_BINARY_DIR}/bin")
  endif()
endif()

if(WIN32) # library to help with demangling variable names on Windows
  target_link_libraries(lupnt PRIVATE Dbghelp)
endif()

install(
  TARGETS lupnt
  EXPORT LPT-exports
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

list(APPEND LPT_EXPORTED_TARGETS lupnt)
set(LPT_EXPORTED_TARGETS "${LPT_EXPORTED_TARGETS}" PARENT_SCOPE)