###############################################################################

# Use bundled Eigen include path.
# Clear any variables set by FindEigen3
if(EIGEN3_INCLUDE_DIR)
    set(EIGEN3_INCLUDE_DIR NOTFOUND CACHE STRING "" FORCE)
endif()

# set full path to be used by external projects
# this will be added to LPT_INCLUDE_DIR by LPT_extra.cmake.in
set(LPT_EIGEN_INCLUDE_FOR_INSTALL "include/lupnt/3rdparty/Eigen/")

# The actual include directory (for BUILD cmake target interface):
set(LPT_EIGEN_INCLUDE_FOR_BUILD "${LPT_SOURCE_DIR}/lupnt/3rdparty/Eigen/")
