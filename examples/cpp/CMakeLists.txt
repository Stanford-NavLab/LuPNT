# define variables
set(groupName "examples")
set(linkLibraries "lupnt")

# Get all script files
file(GLOB script_files "*.cc" "**/*.cc")

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
    message(STATUS "  ${script_name}")
endforeach()

# Create executables
foreach(script_src IN ITEMS ${script_srcs})
    # Get script base name
    get_filename_component(script_name ${script_src} NAME_WE)

    # Add executable
    add_executable(${script_name} ${script_src} ${script_headers})
    target_link_libraries(${script_name} ${linkLibraries})

    # Apply user build flags from CMake cache variables:
    # gtsam_apply_build_flags(${script_name})

    # Add target dependencies
    # add_dependencies(${groupName} ${script_name})

    # Add TOPSRCDIR
    set_property(SOURCE ${script_src} APPEND PROPERTY COMPILE_DEFINITIONS "TOPSRCDIR=\"${LuPNT_SOURCE_DIR}\"")

    # Configure target folder (for MSVC and Xcode)
    set_property(TARGET ${script_name} PROPERTY FOLDER "${groupName}")
endforeach()