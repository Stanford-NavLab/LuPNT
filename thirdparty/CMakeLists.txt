# ---- Add dependencies via CPM ----

# PackageProject
CPMAddPackage(
  NAME PackageProject.cmake
  GITHUB_REPOSITORY TheLartians/PackageProject.cmake
  VERSION 1.11.2
)

# fmt
CPMAddPackage(
  NAME fmt
  GIT_TAG 11.0.2
  GITHUB_REPOSITORY fmtlib/fmt
  OPTIONS "FMT_INSTALL ON" "CMAKE"
)

# Eigen
CPMAddPackage(
  NAME Eigen
  GIT_TAG 3.4
  GITLAB_REPOSITORY libeigen/eigen
  OPTIONS "EIGEN_BUILD_DOC OFF" "EIGEN_BUILD_BTL OFF" "EIGEN_BUILD_PKGCONFIG OFF"
)

# autodiff
CPMAddPackage(
  NAME autodiff
  GIT_TAG v1.1.2
  GITHUB_REPOSITORY autodiff/autodiff
  OPTIONS "AUTODIFF_BUILD_TESTS OFF" "AUTODIFF_BUILD_EXAMPLES OFF" "AUTODIFF_BUILD_DOCS OFF"
          "AUTODIFF_BUILD_PYTHON OFF"
)

# cspice
CPMAddPackage(
  NAME cspice
  GITHUB_REPOSITORY tudat-team/cspice-cmake
  GIT_TAG master
  # DOWNLOAD_ONLY TRUE
)

# include(ExternalProject) # ExternalProject_Add(cspice #     SOURCE_DIR ${cspice_SOURCE_DIR} # )
# ExternalProject_Add(cspice GIT_REPOSITORY https://github.com/tudat-team/cspice-cmake.git GIT_TAG
# master )

# matplotplusplus
CPMAddPackage(
  NAME matplotplusplus
  GIT_TAG v1.2.1
  GITHUB_REPOSITORY alandefreitas/matplotplusplus
  OPTIONS "MATPLOTPP_BUILD_INSTALLER ON" "MATPLOTPP_BUILD_PACKAGE ON"
)
install(
  TARGETS cimg
  EXPORT Matplot++Targets
  INCLUDES
  DESTINATION include
)
install(
  TARGETS nodesoup
  EXPORT Matplot++Targets
  INCLUDES
  DESTINATION include
)

# export cspice_SOURCE_DIR dir to parent scope set(cspice_SOURCE_DIR ${cspice_SOURCE_DIR}
# PARENT_SCOPE) set(matplot_SOURCE_DIR ${matplot_SOURCE_DIR} PARENT_SCOPE)

# --- Create interface for libraries ---
# Note: this disables warnings for third-party code add_library(thirdparty INTERFACE)
# set_target_properties(thirdparty PROPERTIES CXX_STANDARD 17)

# if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
# target_compile_options(autodiff INTERFACE -w) elseif(MSVC) target_compile_options(autodiff
# INTERFACE /W0) endif()

# target_link_libraries(thirdparty INTERFACE SYSTEM fmt::fmt Eigen3::Eigen cspice matplot )
# target_include_directories(thirdparty INTERFACE
# $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>
# $<BUILD_INTERFACE:${cspice_SOURCE_DIR}/include> )

# install(TARGETS thirdparty EXPORT thirdpartyTargets INCLUDES DESTINATION include )

# install(EXPORT thirdpartyTargets FILE thirdpartyTargets.cmake NAMESPACE ${PROJECT_NAME}::
# DESTINATION lib/cmake/${PROJECT_NAME} )

# # Create a variable called LUPNT_DEPENDENCIES that should be available at the parent scope
set(LUPNT_DEPENDENCIES
    "fmt 7.1.3;Eigen 3.4;autodiff 1.1.2;cspice;matplot 1.2.1"
    PARENT_SCOPE
)
