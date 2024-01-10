# License - MIT - Open3D
# Option 1: Do not define "PYTHON_EXECUTABLE", but run `cmake ..` within your
#           virtual environment. CMake will pick up the python executable in the
#           virtual environment.
# Option 2: You can also define `cmake -DPYTHON_EXECUTABLE` to specify a python
#           executable.
set(PYBIND_PACKAGE_FOLDER "${PYPI_PACKAGE_NAME}_pybind")

if (NOT PYTHON_EXECUTABLE)
    # find_program will returns the python executable in current PATH, which
    # works with virtualenv
    find_program(PYTHON_IN_PATH "python")
    set(PYTHON_EXECUTABLE ${PYTHON_IN_PATH})
endif()
message(STATUS "Using Python executable: ${PYTHON_EXECUTABLE}")

# We need to get python version to configure some meta files
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print('%d.%d' % (sys.version_info.major, sys.version_info.minor))"
    OUTPUT_VARIABLE PYTHON_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Using Python version: ${PYTHON_VERSION}")

execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print('%d' % (sys.version_info.major,))"
    OUTPUT_VARIABLE PYTHON_VERSION_MAJOR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Using Python version major: ${PYTHON_VERSION_MAJOR}")

set(PACKAGE_NAME pylupnt_pybind)

file(GLOB_RECURSE PY_ALL_SOURCE_FILES "${PYBIND_PACKAGE_FOLDER}/*.cc")

message(STATUS "PY_ALL_SOURCE_FILES:")
foreach(PY_ALL_SOURCE_FILE ${PY_ALL_SOURCE_FILES})
    message(STATUS "    ${PY_ALL_SOURCE_FILE}")
endforeach()


# Compile pybind cpp files into a Python module
pybind11_add_module(${PACKAGE_NAME}
    ${PY_ALL_SOURCE_FILES}
)

# Custom command to generate Python stubs after building the module
# add_custom_command(TARGET ${PACKAGE_NAME} POST_BUILD
#                    COMMAND ${PYTHON_EXECUTABLE} -m pybind11_stubgen ${PACKAGE_NAME}
#                    COMMENT "Generating Python stubs for ${PACKAGE_NAME}")


# Add the root directory of the project to the include list
target_include_directories(${PACKAGE_NAME} PRIVATE
    ${PROJECT_SOURCE_DIR}
)

# # This is needed to avoid compilation issues when using clang-9 and c++17
if(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 10.0.0)
    target_compile_options(pylupnt PRIVATE -fsized-deallocation)
endif()

# Suppress Pybind11 warnings
target_include_directories(${PACKAGE_NAME} SYSTEM PRIVATE
    ${PYBIND11_INCLUDE_DIR}
)

# Link _pylupnt against lupnt C++ library
target_link_libraries(${PACKAGE_NAME} PRIVATE lupnt)
if (${PYTHON_VERSION_MAJOR} EQUAL 2)
    target_compile_definitions(${PACKAGE_NAME} PRIVATE PYTHON_2_FALLBACK)
endif ()

# At `make`: xxx.so (or the equivalents) will be created at
# PYTHON_COMPILED_MODULE_DIR. The default locaiton is `build/lib/Python`
set(PYTHON_COMPILED_MODULE_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Python")
set_target_properties(${PACKAGE_NAME} PROPERTIES
                      FOLDER "Python"
                      LIBRARY_OUTPUT_DIRECTORY "${PYTHON_COMPILED_MODULE_DIR}"
                      ARCHIVE_OUTPUT_DIRECTORY "${PYTHON_COMPILED_MODULE_DIR}")


                      
# Use `make python-package` to create the python package in the build directory
# The python package will be created at PYTHON_PACKAGE_DIR. It contains:
# 1) Pure-python code and misc files, copied from src/Python/package
# 2) The compiled python-C++ module, i.e. xxx.so (or the equivalents)
# 3) Configured files and supporting files
# Note: `make python-package` clears PYTHON_COMPILED_MODULE_DIR first every time
set(PYTHON_PACKAGE_DST_DIR "${CMAKE_BINARY_DIR}/lib/python_package")
message(STATUS "PYPI_PACKAGE_NAME: ${PYPI_PACKAGE_NAME}")

macro(print_all_variables)
    message(STATUS "print_all_variables------------------------------------------{")
    get_cmake_property(_variableNames VARIABLES)
    foreach (_variableName ${_variableNames})
        message(STATUS "${_variableName}=${${_variableName}}")
    endforeach()
    message(STATUS "print_all_variables------------------------------------------}")
endmacro()

# print_all_variables()

add_custom_target(python-package
    COMMAND ${CMAKE_COMMAND}
            -DPYTHON_PACKAGE_SRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DPYTHON_PACKAGE_DST_DIR=${PYTHON_PACKAGE_DST_DIR}
            -DPYTHON_VERSION=${PYTHON_VERSION}
            -DPYTHON_COMPILED_MODULE_PATH=$<TARGET_FILE:${PACKAGE_NAME}>
            -DORIG_CMAKE_LIBRARY_OUTPUT_DIRECTORY=${ORIG_CMAKE_LIBRARY_OUTPUT_DIRECTORY}
            -DPROJECT_EMAIL=${PROJECT_EMAIL}
            -DPROJECT_HOME=${PROJECT_HOME}
            -DPROJECT_DOCS=${PROJECT_DOCS}
            -DPROJECT_CODE=${PROJECT_CODE}
            -DPROJECT_ISSUES=${PROJECT_ISSUES}
            -DPROJECT_VERSION=${PROJECT_VERSION}
            -DPYPI_PACKAGE_NAME=${PYPI_PACKAGE_NAME}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/make_python_package.cmake
)

# Use `make pip-package` to create the pip package in the build directory
add_custom_target(pip-package
    COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_wheel --dist-dir pip_package
    COMMAND echo "pip wheel created at ${PYTHON_PACKAGE_DST_DIR}/pip_package"
    WORKING_DIRECTORY ${PYTHON_PACKAGE_DST_DIR}
    DEPENDS python-package
)

# Use `make install-pip-package` to install pip wheel package to the current
# python environment.
add_custom_target(install-pip-package
    COMMAND ${CMAKE_COMMAND}
            -DPYTHON_PACKAGE_DST_DIR=${PYTHON_PACKAGE_DST_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/make_install_pip_package.cmake
    DEPENDS pip-package
)





# # define variables
# set(groupName "bindings")

# # Collect all source files from the current directory for the Python module
# file(GLOB_RECURSE CXX_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.cc")

# # Compile pybind cpp files into a Python module
# pybind11_add_module(_pylupnt ${CXX_FILES})

# # Add the root directory of the project to the include list
# target_include_directories(_pylupnt PRIVATE ${PROJECT_SOURCE_DIR})

# # Link _pylupnt against lupnt C++ library
# target_link_libraries(_pylupnt PRIVATE lupnt)

# # This is needed to avoid compilation issues when using clang-9 and c++17
# if(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 10.0.0)
#     target_compile_options(_pylupnt PRIVATE -fsized-deallocation)
# endif()

