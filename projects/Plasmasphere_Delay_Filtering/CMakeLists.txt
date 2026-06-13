cmake_minimum_required(VERSION 3.15)
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON)

# Set project name as the name of the current source directory
get_filename_component(project_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${project_name} LANGUAGES CXX)

file(
  GLOB main_files
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  main*.cc
)

# Collect all source files from src/ subdirectory
file(GLOB_RECURSE src_files ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

# ---- Dependencies ----
include(${CMAKE_CURRENT_LIST_DIR}/../../cmake/CPM.cmake)

# Add TBB ------ Make conda discoverable for find_package
list(PREPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
set(TBB_DIR "$ENV{CONDA_PREFIX}/lib/cmake/TBB")
find_package(TBB CONFIG REQUIRED COMPONENTS tbb tbbmalloc)

foreach(main_file ${main_files})
  # Keep the name with the main_ prefix, i.e., use main_dataset for target
  get_filename_component(target_name ${main_file} NAME_WE)

  add_executable(${target_name} ${main_file} ${src_files})
  target_include_directories(${target_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(${target_name} PRIVATE lupnt::lupnt TBB::tbb TBB::tbbmalloc)

  # Ensure the linker searches your env’s lib dir (fixes -ltbb)
  target_link_directories(${target_name} PRIVATE "$ENV{CONDA_PREFIX}/lib")

  # if MACOS
  if(APPLE)
    set_target_properties(
      ${target_name}
      PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE
                 INSTALL_RPATH "@loader_path;@loader_path/../lib;$ENV{CONDA_PREFIX}/lib"
    )
  endif()

  set_target_properties(${target_name} PROPERTIES CXX_STANDARD 20)
  set_target_properties(
    ${target_name}
    PROPERTIES FOLDER "projects/${project_name}"
               RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/projects/${project_name}
  )

  # Path to the built shared library
  set(LUPNT_LIB "${CMAKE_BINARY_DIR}/lupnt/liblupnt.dylib")

  # Copy liblupnt.dylib into the executable directory after build
  add_custom_command(
    TARGET ${target_name}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${LUPNT_LIB}"
            "$<TARGET_FILE_DIR:${target_name}>/"
    COMMENT "Copying liblupnt.dylib to output directory"
  )
endforeach()
