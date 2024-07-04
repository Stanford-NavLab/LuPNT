# Configure the setup.py file
configure_file(setup.py.in ${CMAKE_CURRENT_BINARY_DIR}/setup.py)

# Create a custom target to build the python package during build stage
add_custom_target(lupnt-setuptools ALL
    DEPENDS pylupnt
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/lupnt ${CMAKE_CURRENT_BINARY_DIR}/lupnt
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pylupnt> ${CMAKE_CURRENT_BINARY_DIR}/lupnt
    COMMAND ${PYTHON_EXECUTABLE} setup.py --quiet build --force
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/lupnt/$<TARGET_FILE_NAME:pylupnt>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Ensure the path where the python package is installed is not empty
if(NOT DEFINED LUPNT_PYTHON_INSTALL_PREFIX)
    file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX} LUPNT_PYTHON_INSTALL_PREFIX)
endif()

# Create an install target for the lupnt python package
install(CODE "execute_process(COMMAND pip install . --no-deps --force-reinstall --no-compile --prefix=${LUPNT_PYTHON_INSTALL_PREFIX} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")
