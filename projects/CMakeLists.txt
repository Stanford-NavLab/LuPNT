cmake_minimum_required(VERSION 3.15)
project(MainProject)

# Get all immediate subdirectories
file(
  GLOB subdirs
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*
)

# Iterate through the subdirectories
foreach(subdir ${subdirs})
  if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${subdir})
    # Check if the subdirectory contains a CMakeLists.txt file
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt)
      # Add the subdirectory to the build
      add_subdirectory(${subdir})
      message(STATUS "Added subdirectory: ${subdir}")
    endif()
  endif()
endforeach()
