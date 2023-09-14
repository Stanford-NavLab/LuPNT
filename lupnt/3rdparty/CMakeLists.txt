# install Eigen -------------------------------------------------------
# # Find plain .h files
# file(GLOB_RECURSE eigen_headers "${CMAKE_CURRENT_SOURCE_DIR}/Eigen/Eigen/*.h")

# # Find header files without extension
# file(GLOB eigen_dir_headers_all "Eigen/Eigen/*")
# foreach(eigen_dir ${eigen_dir_headers_all})
#     get_filename_component(filename ${eigen_dir} NAME)
#     if (NOT ((${filename} MATCHES "CMakeLists.txt") OR (${filename} MATCHES "src") OR (${filename} MATCHES ".svn")))
#         list(APPEND eigen_headers "${CMAKE_CURRENT_SOURCE_DIR}/Eigen/Eigen/${filename}")
#         install(FILES Eigen/Eigen/${filename} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lupnt/3rdparty/Eigen/Eigen)
#     endif()
# endforeach(eigen_dir)

# # Add to project source
# set(eigen_headers ${eigen_headers} PARENT_SCOPE)

# # install Eigen - only the headers in our 3rdparty directory
# install(DIRECTORY Eigen/Eigen
#     DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lupnt/3rdparty/Eigen
#     FILES_MATCHING PATTERN "*.h")

# install autodiff -----------------------------------------------------
# Find plain .h files
file(GLOB_RECURSE autodiff_headers "${CMAKE_CURRENT_SOURCE_DIR}/autodiff/autodiff/*.hpp")

# Add to project source
set(autodiff_headers ${autodiff_headers} PARENT_SCOPE)

# install autodiff - only the headers in our 3rdparty directory
install(DIRECTORY autodiff/autodiff
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lupnt/3rdparty/autodiff
    FILES_MATCHING PATTERN "*.hpp")

# install cspice ------------------------------------------------------------
file(GLOB spice_headers "{CMAKE_CURRENT_SOURCE_DIR}/cspice/*h")

# Add to project source
set(spice_headers ${spice_headers} PARENT_SCOPE)

# install cspice - only the headers in our 3rdparty directory
install(DIRECTORY cspice
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lupnt/3rdparty/cspice
    FILES_MATCHING PATTERN "*.h")