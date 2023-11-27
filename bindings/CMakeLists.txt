# define variables
set(groupName "bindings")

# Collect all source files from the current directory for the Python module
file(GLOB_RECURSE CXX_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.cc")

# Compile pybind cpp files into a Python module
pybind11_add_module(pylupnt ${CXX_FILES})

# Add the root directory of the project to the include list
target_include_directories(pylupnt PRIVATE ${PROJECT_SOURCE_DIR})

# Link pylupnt against pylupnt C++ library
target_link_libraries(pylupnt PRIVATE lupnt)

# This is needed to avoid compilation issues when using clang-9 and c++17
if(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 10.0.0)
    target_compile_options(pylupnt PRIVATE -fsized-deallocation)
endif()

# copy to PYTHON_LIBRARIES
add_custom_command(TARGET pylupnt POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pylupnt> ${PYTHON_LIBRARIES}
        COMMENT "Copying pylupnt to ${PYTHON_LIBRARIES}"
        VERBATIM
)