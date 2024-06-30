project(lupnt LANGUAGES CXX)

# We split the library in to separate subfolders, each containing
# tests, timing, and an optional convenience library.
# The following variable is the master list of subdirs to add
set(LuPNT_subdirs
  core
  datasets
  dynamics
  numerics
  agents
  physics
  measurements
)

set(LuPNT_srcs)

# Add the "lupnt" directory to the include search paths
include_directories("${CMAKE_SOURCE_DIR}/lupnt")

# assemble core libaries
message(STATUS "Configuring modules")
foreach(subdir ${LuPNT_subdirs})
  # Build convenience libraries

  # Build local library and tests
  message(STATUS "   ${subdir}")

  # Include header files so they show up in Visual Studio
  file(GLOB_RECURSE subdir_srcs "${subdir}/*.cc" "${subdir}/*.h") 
  # Print the source files (only the filename, not the full path)
  foreach(file ${subdir_srcs})
    get_filename_component(file_name ${file} NAME)
    message(STATUS "     ${file_name}")
  endforeach(file)
  
  set(${subdir}_srcs ${subdir_srcs})

  add_subdirectory(${subdir})
endforeach(subdir)

# To add additional sources to gtsam when building the full library (static or shared)
# append the subfolder with _srcs appended to the end to this list
set(LuPNT_srcs ${third_party_srcs})

foreach(subdir ${LuPNT_subdirs})
  list(APPEND LuPNT_srcs ${${subdir}_srcs})
endforeach(subdir)

# build shared and static versions of the library
message(STATUS "Configuring lupnt - shared: ${BUILD_SHARED_LIBS}")

# BUILD_SHARED_LIBS automatically defines static/shared libs:
add_library(lupnt ${LuPNT_srcs})

# add cspice
target_link_libraries(lupnt PUBLIC "${CMAKE_SOURCE_DIR}/third_party/cspice/cspice.a")

# add MatPlot++
target_link_libraries(lupnt PUBLIC matplot)

# add OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(lupnt PUBLIC OpenMP::OpenMP_CXX)

# add Boost
target_link_libraries(lupnt PUBLIC ${Boost_LIBRARIES})

# include libraries
target_include_directories(lupnt BEFORE PUBLIC # main gtsam includes:
  $<BUILD_INTERFACE:${LuPNT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/>
)

# third_party libraries: use the "system" flag so they are included via "-isystem"
# and warnings (and warnings-considered-errors) in those headers are not
# reported as warnings/errors in our targets:
target_include_directories(lupnt SYSTEM BEFORE PUBLIC

  # autodiff
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/third_party/autodiff>
  $<INSTALL_INTERFACE:include/lupnt/third_party/autodiff>

  # Eigen
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/third_party/eigen>
  $<INSTALL_INTERFACE:${include/lupnt/third_party/eigen}>

  # cspice
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/third_party/cspice>
  $<INSTALL_INTERFACE:include/lupnt/third_party/cspice>
)

if(WIN32) # Add 'lib' prefix to static library to avoid filename collision with shared library
  if(NOT BUILD_SHARED_LIBS)
    set_target_properties(lupnt PROPERTIES
      PREFIX "lib"
      COMPILE_DEFINITIONS LuPNT_IMPORT_STATIC)
  else()
    set_target_properties(lupnt PROPERTIES
      PREFIX ""
      DEFINE_SYMBOL LuPNT_EXPORTS
      RUNTIME_OUTPUT_DIRECTORY "${LuPNT_BINARY_DIR}/bin")
  endif()
endif()

if(WIN32) # library to help with demangling variable names on Windows
  target_link_libraries(lupnt PRIVATE Dbghelp)
endif()

install(
  TARGETS lupnt
  EXPORT LuPNT-exports
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

list(APPEND LuPNT_EXPORTED_TARGETS lupnt)
set(LuPNT_EXPORTED_TARGETS "${LuPNT_EXPORTED_TARGETS}" PARENT_SCOPE)


# *****************************************************************************
# Create a single header file lupnt.h with all includes

# Navigate to the lupnt directory
set(LUPNT_DIR "${CMAKE_SOURCE_DIR}/lupnt")

# Initialize or clear the lupnt.h file content before appending
file(WRITE "${LUPNT_DIR}/lupnt.h" "")

# Function to process directories and create lupnt.h
function(process_directory dir)
    # Find all subdirectories
    file(GLOB SUBDIRS LIST_DIRECTORIES true RELATIVE ${dir} ${dir}/*)
    list(FILTER SUBDIRS EXCLUDE REGEX "^\\.|\\.\\./")

    foreach(subdir ${SUBDIRS})
        if(IS_DIRECTORY "${dir}/${subdir}")
            # Format the directory name for a comment
            string(REGEX REPLACE "^${dir}/" "" subdir_name ${subdir})
            file(APPEND "${LUPNT_DIR}/lupnt.h" "// ${subdir_name}\n")

            # Find all .h files in the current subdirectory
            file(GLOB HEADER_FILES RELATIVE ${dir}/${subdir} ${dir}/${subdir}/*.h)
            list(SORT HEADER_FILES)

            foreach(file ${HEADER_FILES})
                # Create relative include path
                file(APPEND "${LUPNT_DIR}/lupnt.h" "#include \"${subdir}/${file}\"\n")
            endforeach()

            # Add an empty line for readability between subdirectory sections
            file(APPEND "${LUPNT_DIR}/lupnt.h" "\n")
        endif()
    endforeach()
endfunction()

# Process the lupnt directory
process_directory(${LUPNT_DIR})

message(STATUS "A single header file lupnt.h has been created with all includes.")