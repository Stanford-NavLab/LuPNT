# This cmake module determines if a conda environment is active.
#
# Note: CMAKE_INSTALL_PREFIX is set to the environment variable CONDA_PREFIX if a conda environment
# is found active.
#
# This automatic behavior can be overriden by manually specifying a different CMAKE_INSTALL_PREFIX.

# Skip if this file has already been included
if(CONDA_AWARE_INCLUDED)
  return()
else()
  set(CONDA_AWARE_INCLUDED TRUE)
endif()

# Check if a conda environment is active
if(DEFINED ENV{CONDA_PREFIX})
  # Show that conda env has been recognized
  message(STATUS "CondaAware: Conda environment detected!")
  message(STATUS "CondaAware: Found environment variable CONDA_PREFIX=$ENV{CONDA_PREFIX}")

  # Check if environment variable PYTHON is defined, and if so, set PYTHON_EXECUTABLE to PYTHON
  if(DEFINED ENV{PYTHON})
    message(STATUS "CondaAware: Found environment variable PYTHON=$ENV{PYTHON}")
    message(STATUS "CondaAware: Setting PYTHON_EXECUTABLE to PYTHON=$ENV{PYTHON}")
    set(PYTHON_EXECUTABLE $ENV{PYTHON})
  endif()

  # Ensure python executable is properly selected if not specified yet using cmake variables
  # PYTHON_EXECUTABLE or environment variable PYTHON. The code below tries to set PYTHON_EXECUTABLE
  # with the python executable in the activated conda environment and not in the base environment
  # given by the environment variable CONDA_PYTHON_EXE!
  if(NOT DEFINED PYTHON_EXECUTABLE)
    if(UNIX)
      message(STATUS "CondaAware: Setting PYTHON_EXECUTABLE=$ENV{CONDA_PREFIX}/bin/python")
      set(PYTHON_EXECUTABLE "$ENV{CONDA_PREFIX}/bin/python")
    endif()

    if(WIN32)
      message(STATUS "CondaAware: Setting PYTHON_EXECUTABLE=$ENV{CONDA_PREFIX}\\python.exe")
      set(PYTHON_EXECUTABLE "$ENV{CONDA_PREFIX}\\python.exe")
    endif()

    if(NOT DEFINED PYTHON_EXECUTABLE)
      message(
        FATAL_ERROR
          "CondaAware: Could not determine a value for PYTHON_EXECUTABLE. Expecting Unix or Windows systems."
      )
    endif()
  endif()

  # Set auxiliary variable CONDA_AWARE_PREFIX as according to the logic below.
  if(DEFINED ENV{CONDA_BUILD})
    message(STATUS "CondaAware: Detected conda build task (e.g., in a conda-forge build)!")
    if(UNIX)
      set(CONDA_AWARE_PREFIX "$ENV{PREFIX}")
    endif()
    if(WIN32)
      set(CONDA_AWARE_PREFIX "$ENV{LIBRARY_PREFIX}")
    endif()
  else()
    if(UNIX)
      set(CONDA_AWARE_PREFIX "$ENV{CONDA_PREFIX}")
    endif()

    if(WIN32)
      set(CONDA_AWARE_PREFIX "$ENV{CONDA_PREFIX}\\Library")
    endif()
  endif()

  # Check if CONDA_AWARE_PREFIX has been successfully set
  if(DEFINED CONDA_AWARE_PREFIX)
    message(STATUS "CondaAware: Set CONDA_AWARE_PREFIX=${CONDA_AWARE_PREFIX}")
  else()
    message(
      FATAL_ERROR
        "CondaAware: Could not determine a value for CONDA_AWARE_PREFIX. Expecting Unix or Windows systems."
    )
  endif()

  # Set CMAKE_INSTALL_PREFIX to CONDA_AWARE_PREFIX if not specified by the user
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(
      STATUS "CondaAware: Setting CMAKE_INSTALL_PREFIX=CONDA_AWARE_PREFIX=${CONDA_AWARE_PREFIX}"
    )
    set(CMAKE_INSTALL_PREFIX ${CONDA_AWARE_PREFIX})
  endif()

  # Ensure dependencies from the conda environment are used instead of those from the system.
  # Prioritize conda libraries by putting them first in CMAKE_PREFIX_PATH
  list(PREPEND CMAKE_PREFIX_PATH ${CONDA_AWARE_PREFIX})
  message(STATUS "CondaAware: Prepended ${CONDA_AWARE_PREFIX} to CMAKE_PREFIX_PATH")

  # Ensure include directory in conda environment is known to the project. Wrapped in
  # $<BUILD_INTERFACE:> so it is never baked into install(EXPORT ...) interfaces for targets whose
  # CMAKE_INSTALL_PREFIX-relative source path (e.g. an in-tree .pixi/envs/.../include) would
  # otherwise trip CMake's "path is prefixed in the source directory" export check.
  include_directories("$<BUILD_INTERFACE:${CONDA_AWARE_PREFIX}/include>")
  message(STATUS "CondaAware: Appended ${CONDA_AWARE_PREFIX}/include to include directories")

  # Ensure library directory in conda environment is known to the project
  link_directories(${CONDA_AWARE_PREFIX}/lib)
  message(STATUS "CondaAware: Appended ${CONDA_AWARE_PREFIX}/lib to link directories")

  # Set RPATH to ensure conda environment libraries are used at runtime Skip for Python modules as
  # pybind11 handles RPATH automatically
  if(NOT DEFINED CMAKE_CURRENT_SOURCE_DIR OR NOT CMAKE_CURRENT_SOURCE_DIR MATCHES
                                             ".*python.*bindings.*"
  )
    # Check if our path is already in the RPATH (either as a single path or in a list)
    set(rpath_already_set FALSE)

    # Check CMAKE_INSTALL_RPATH
    if(DEFINED CMAKE_INSTALL_RPATH)
      if(CMAKE_INSTALL_RPATH STREQUAL "${CONDA_AWARE_PREFIX}/lib")
        set(rpath_already_set TRUE)
      else()
        # Check if it's in a list
        list(FIND CMAKE_INSTALL_RPATH "${CONDA_AWARE_PREFIX}/lib" rpath_index)
        if(NOT rpath_index EQUAL -1)
          set(rpath_already_set TRUE)
        endif()
      endif()
    endif()

    # Check CMAKE_BUILD_RPATH
    if(NOT rpath_already_set AND DEFINED CMAKE_BUILD_RPATH)
      if(CMAKE_BUILD_RPATH STREQUAL "${CONDA_AWARE_PREFIX}/lib")
        set(rpath_already_set TRUE)
      else()
        # Check if it's in a list
        list(FIND CMAKE_BUILD_RPATH "${CONDA_AWARE_PREFIX}/lib" rpath_index)
        if(NOT rpath_index EQUAL -1)
          set(rpath_already_set TRUE)
        endif()
      endif()
    endif()

    if(rpath_already_set)
      message(
        STATUS
          "CondaAware: ${CONDA_AWARE_PREFIX}/lib already in RPATH (CMAKE_INSTALL_RPATH or CMAKE_BUILD_RPATH)"
      )
    else()
      if(NOT DEFINED CMAKE_INSTALL_RPATH OR CMAKE_INSTALL_RPATH STREQUAL "")
        set(CMAKE_INSTALL_RPATH "${CONDA_AWARE_PREFIX}/lib")
        message(STATUS "CondaAware: Set CMAKE_INSTALL_RPATH=${CONDA_AWARE_PREFIX}/lib")
      else()
        list(APPEND CMAKE_INSTALL_RPATH "${CONDA_AWARE_PREFIX}/lib")
        message(
          STATUS "CondaAware: Appended ${CONDA_AWARE_PREFIX}/lib to existing CMAKE_INSTALL_RPATH"
        )
      endif()
    endif()
  else()
    message(STATUS "CondaAware: Skipping RPATH setup for Python module build")
  endif()

  # For macOS and when using Ninja generator, also set the build RPATH
  if(APPLE OR CMAKE_GENERATOR STREQUAL "Ninja")
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif()

  set(CMAKE_FIND_ROOT_PATH ${CONDA_AWARE_PREFIX})
  set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

endif()
