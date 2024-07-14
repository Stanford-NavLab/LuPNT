# Create the lupnt static library
file(GLOB_RECURSE LUPNT_SOURCES CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/lupnt/*.cc")
file(GLOB_RECURSE LUPNT_HEADERS CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/lupnt/*.h")

# Create a single header file lupnt.h with all includes
set(LUPNT_DIR "${CMAKE_SOURCE_DIR}/lupnt")
file(WRITE "${LUPNT_DIR}/lupnt.h" "")
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
process_directory(${LUPNT_DIR})

# Add an alias lupnt::lupnt to the target library lupnt
add_library(lupnt STATIC ${LUPNT_SOURCES} ${LUPNT_HEADERS})
add_library(lupnt::lupnt ALIAS lupnt)
# Set lupnt compilation features to be propagated to client code.
target_compile_features(lupnt PUBLIC cxx_std_23)

# Find CSPICE include directory and library
find_path(CSPICE_INCLUDE_DIR NAMES SpiceUsr.h
          PATHS
          /opt/homebrew/opt/cspice/include
          /usr/local/include
          /usr/include
          $ENV{CSPICE_HOME}/include
          C:/cspice/include
          DOC "Path to the CSPICE include directory")
find_library(CSPICE_LIBRARY NAMES cspice
             PATHS
             /opt/homebrew/opt/cspice/lib
             /usr/local/lib
             /usr/lib
             $ENV{CSPICE_HOME}/lib
             C:/cspice/lib
             DOC "Path to the CSPICE library")
if (NOT CSPICE_INCLUDE_DIR OR NOT CSPICE_LIBRARY)
    message(FATAL_ERROR "CSPICE library not found. Please install CSPICE.")
endif()

# Path to the OpenMP include directory
find_path(OpenMP_INCLUDE_DIR NAMES omp.h
          PATHS
          /opt/homebrew/opt/libomp/include
          /usr/local/opt/libomp/include
          /usr/local/include
          /usr/include
          $ENV{OpenMP_ROOT}/include
          DOC "Path to the OpenMP include directory")
find_library(OpenMP_omp_LIBRARY NAMES omp
             PATHS
             /opt/homebrew/opt/libomp/lib
             /usr/local/opt/libomp/lib
             $ENV{OpenMP_ROOT}/lib
             DOC "Path to the CSPICE library")
if (NOT OpenMP_INCLUDE_DIR OR NOT OpenMP_omp_LIBRARY)
    message(FATAL_ERROR "OpenMP library not found. Please install OpenMP.")
endif()

target_include_directories(lupnt
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        ${CSPICE_INCLUDE_DIR} # Add CSPICE include directory
        ${OpenMP_INCLUDE_DIR} # Add OpenMP include directory
)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(autodiff 1 REQUIRED NO_MODULE)
find_package(Matplot++ 1.2 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)
target_link_libraries(lupnt PUBLIC 
    Eigen3::Eigen
    autodiff::autodiff
    Matplot++::matplot
    OpenMP::OpenMP_CXX
    ${CSPICE_LIBRARY} # Link against CSPICE library
)

# Install lupnt static library
install(TARGETS lupnt EXPORT lupntTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install lupnt header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/lupnt
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers
    PATTERN "CMakeLists.txt" EXCLUDE)