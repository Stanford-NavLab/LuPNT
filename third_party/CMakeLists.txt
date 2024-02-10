# Autodiff --------------------------------------------------------------------

# Find plain .h files
file(GLOB_RECURSE autodiff_headers "${CMAKE_CURRENT_SOURCE_DIR}/autodiff/autodiff/*.hpp")

# Add to project source
set(autodiff_headers ${autodiff_headers} PARENT_SCOPE)

# install autodiff - only the headers in our third_party directory
install(DIRECTORY autodiff/autodiff
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lupnt/third_party/autodiff
    FILES_MATCHING PATTERN "*.hpp")

# Assert headers are found
if(NOT autodiff_headers)
    message(FATAL_ERROR "Autodiff headers not found")
endif()

# cspice ----------------------------------------------------------------------

file(GLOB_RECURSE cspice_headers "${CMAKE_CURRENT_SOURCE_DIR}/cspice/include/*.h")

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cspice/cspice")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/cspice/cspice")
    foreach(header ${cspice_headers})
        file(COPY ${header} DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/cspice/cspice")
    endforeach()
endif()


file(GLOB_RECURSE cspice_headers "${CMAKE_CURRENT_SOURCE_DIR}/cspice/cspice/*.h")

# Add to project source
set(cspice_headers ${cspice_headers} PARENT_SCOPE)

# install cspice - only the headers in our third_party directory
install(DIRECTORY cspice/include
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lupnt/third_party/cspice
    FILES_MATCHING PATTERN "*.h")

if(NOT cspice_headers)
    message(FATAL_ERROR "cspice headers not found")
endif()

# Boost -----------------------------------------------------------------------

# Store these in variables so they are automatically replicated in GTSAMConfig.cmake and such.
set(BOOST_FIND_MINIMUM_VERSION 1.65)
set(BOOST_FIND_MINIMUM_COMPONENTS serialization system filesystem thread program_options date_time timer chrono regex)

find_package(Boost ${BOOST_FIND_MINIMUM_VERSION} COMPONENTS ${BOOST_FIND_MINIMUM_COMPONENTS})

# Required components
if(NOT Boost_SERIALIZATION_LIBRARY OR NOT Boost_SYSTEM_LIBRARY OR NOT Boost_FILESYSTEM_LIBRARY OR
    NOT Boost_THREAD_LIBRARY OR NOT Boost_DATE_TIME_LIBRARY)
  message(FATAL_ERROR "Missing required Boost components >= v1.65, please install/upgrade Boost or configure your search paths.")
else()
    message("Boost lib found")
endif()

# Python ----------------------------------------------------------------------
# set(PYBIND11_FINDPYTHON ON CACHE BOOL "Force pybind11 to use FindPython") 
add_subdirectory(pybind11)


# Matplot++ -------------------------------------------------------------------
add_subdirectory(matplotplusplus)

# AI-Toolbox -------------------------------------------------------------------

set(LPSOLVE_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen")
add_subdirectory(AI-Toolbox)
file(GLOB_RECURSE autodiff_headers "${CMAKE_CURRENT_SOURCE_DIR}/autodiff/autodiff/*.hpp")

# ------------------------------------------------------------

set(third_party_srcs
  ${eigen_headers} 
  ${autodiff_headers}
  ${cspice_headers}
)