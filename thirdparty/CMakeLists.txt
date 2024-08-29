# ---- Add dependencies via CPM ----

# PackageProject
CPMAddPackage(
  NAME PackageProject.cmake
  GITHUB_REPOSITORY TheLartians/PackageProject.cmake
  VERSION 1.11.2
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
  GITHUB_REPOSITORY gcasadesus/cspice
  GIT_TAG master
  OPTIONS "CMAKE_BUILD_TYPE Release"
)

# matplotplusplus
CPMAddPackage(
  NAME matplotplusplus
  GIT_TAG v1.2.1
  GITHUB_REPOSITORY alandefreitas/matplotplusplus
  OPTIONS "MATPLOTPP_BUILD_INSTALLER ON" "MATPLOTPP_BUILD_PACKAGE ON" "CMAKE_BUILD_TYPE Release"
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

# HighFive
CPMAddPackage(
  NAME HighFive
  GIT_TAG v2.10.0
  GITHUB_REPOSITORY BlueBrain/HighFive
  OPTIONS "HIGHFIVE_UNIT_TESTS OFF" "HIGHFIVE_EXAMPLES OFF" "HIGHFIVE_BUILD_DOCS OFF"
)

set_target_properties(cspice PROPERTIES FOLDER thirdparty)
set_target_properties(matplot PROPERTIES FOLDER thirdparty)
set_target_properties(nodesoup PROPERTIES FOLDER thirdparty)
set_target_properties(demos PROPERTIES FOLDER thirdparty)
set_target_properties(uninstall PROPERTIES FOLDER thirdparty)

# LuPNT dependencies used for installing
set(LUPNT_DEPENDENCIES
    "Eigen3 3.4;autodiff 1.1;cspice;Matplot++ 1.2;OpenMP;HighFive 2.10;Boost"
    PARENT_SCOPE
)
