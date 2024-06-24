# Define variables
set(groupName "examples")
# Enable OpenMP
find_package(OpenMP REQUIRED)
set(linkLibraries "lupnt" OpenMP::OpenMP_CXX)
# Get all script files
file(GLOB_RECURSE script_files "example_*.cc" "**/example_*.cc")

# Separate into source files and headers (allows for adding headers to show up in
# MSVC and Xcode projects).
set(script_srcs "")
set(script_headers "")
foreach(script_file IN ITEMS ${script_files})
    get_filename_component(script_ext ${script_file} EXT)
    if(script_ext MATCHES "(h|H)")
        list(APPEND script_headers ${script_file})
    else()
        list(APPEND script_srcs ${script_file})
    endif()
endforeach()

message(STATUS "Configuring examples")

# Don't put test files in folders in MSVC and Xcode because they're already grouped
source_group("" FILES ${script_srcs} ${script_headers})

# Print script_srcs (only filenames, not full paths)
foreach(script_src IN ITEMS ${script_srcs})
    get_filename_component(script_name ${script_src} NAME)
    message(STATUS "    ${script_name}")
endforeach()

# Create a map to hold subdirectory-based groups
set(subdirectory_groups)

# Create executables and group by subdirectory
foreach(script_src IN ITEMS ${script_srcs})
    # Get script base name and relative subdirectory
    get_filename_component(script_name ${script_src} NAME_WE)
    file(RELATIVE_PATH script_dir ${CMAKE_CURRENT_SOURCE_DIR} ${script_src})
    get_filename_component(script_dir ${script_dir} DIRECTORY)

    # Create a target name as the filename
    set(script_target ${script_name})

    # Add executable
    add_executable(${script_target} ${script_src} ${script_headers})
    target_link_libraries(${script_target} ${linkLibraries})

    # Apply user build flags from CMake cache variables:
    # gtsam_apply_build_flags(${script_target})

    # Add target dependencies
    # add_dependencies(${groupName} ${script_target})

    # Add TOPSRCDIR
    set_property(SOURCE ${script_src} APPEND PROPERTY COMPILE_DEFINITIONS "TOPSRCDIR=\"${LuPNT_SOURCE_DIR}\"")

    # Configure target folder (for MSVC and Xcode)
    set_property(TARGET ${script_target} PROPERTY FOLDER "${groupName}/${script_dir}")

    # Add to subdirectory group
    list(APPEND subdirectory_groups_${script_dir} ${script_target})
endforeach()

# Optionally, print subdirectory groups for debugging
foreach(subdir IN LISTS subdirectory_groups)
    message(STATUS "Group: ${subdir}")
    foreach(target IN LISTS ${subdir})
        message(STATUS "    ${target}")
    endforeach()
endforeach()