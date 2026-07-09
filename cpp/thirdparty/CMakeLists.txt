# ---- Add dependencies via CPM ----

if(NOT DEFINED CPM_SOURCE_CACHE)
  set(CPM_SOURCE_CACHE "$ENV{HOME}/.cache/CPM")
  message(STATUS "Setting CPM_SOURCE_CACHE to ${CPM_SOURCE_CACHE}")
  file(MAKE_DIRECTORY ${CPM_SOURCE_CACHE})
else()
  message(STATUS "CPM_SOURCE_CACHE is already set to ${CPM_SOURCE_CACHE}")
endif()

# PackageProject
CPMAddPackage(
  NAME PackageProject.cmake
  GITHUB_REPOSITORY TheLartians/PackageProject.cmake
  VERSION 1.11.2
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

# matplotplusplus Note: nvpl warnings from FindBLAS/FindLAPACK are harmless - nvpl is only for
# NVIDIA Grace ARM CPUs, not x86_64 These warnings occur inside matplotplusplus's CMakeLists and
# cannot be easily suppressed from here. The warnings are informational and don't affect
# functionality (OpenBLAS is used instead).
CPMAddPackage(
  NAME matplotplusplus
  GIT_TAG v1.2.2
  GITHUB_REPOSITORY alandefreitas/matplotplusplus
  OPTIONS "MATPLOTPP_BUILD_INSTALLER ON" "MATPLOTPP_BUILD_PACKAGE ON" "CMAKE_BUILD_TYPE Release"
)

# Fix macOS SDK compatibility issues with conda clang for matplot
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND APPLE)
  target_compile_definitions(matplot PRIVATE _LIBCPP_DISABLE_AVAILABILITY)
  target_compile_options(
    matplot PRIVATE -Wno-error=availability -Wno-error=deprecated-non-prototype
  )
endif()

# cmake 4.x rejects bare source-tree paths in INTERFACE_INCLUDE_DIRECTORIES of exported targets
# (hard error, not a warning).  When the pixi/conda environment lives inside the project source tree
# (.pixi/), matplotplusplus's 3rd_party cmake adds JPEG/TIFF/PNG/ZLIB/LAPACK/BLAS include dirs to
# the cimg INTERFACE target without $<BUILD_INTERFACE:...> guards.  Wrap any bare paths here so the
# generate-time export check only sees generator expressions.
if(TARGET cimg)
  get_target_property(_cimg_incs cimg INTERFACE_INCLUDE_DIRECTORIES)
  if(_cimg_incs)
    set(_cimg_fixed "")
    foreach(_inc IN LISTS _cimg_incs)
      if(_inc MATCHES "^\\$<")
        list(APPEND _cimg_fixed "${_inc}")
      else()
        list(APPEND _cimg_fixed "$<BUILD_INTERFACE:${_inc}>")
      endif()
    endforeach()
    set_target_properties(cimg PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${_cimg_fixed}")
  endif()
endif()

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
  GIT_TAG v3.0.0-beta2
  GITHUB_REPOSITORY highfive-devs/HighFive
  OPTIONS "HIGHFIVE_UNIT_TESTS OFF" "HIGHFIVE_EXAMPLES OFF" "HIGHFIVE_BUILD_DOCS OFF"
)
target_link_libraries(${PROJECT_NAME} PUBLIC HighFive)

# magic_enum
CPMAddPackage(
  NAME magic_enum
  GIT_TAG v0.9.7
  GITHUB_REPOSITORY Neargye/magic_enum
  OPTIONS "MAGIC_ENUM_OPT_BUILD_EXAMPLES OFF" "MAGIC_ENUM_OPT_BUILD_TESTS OFF"
          "MAGIC_ENUM_OPT_INSTALL ON"
)
install(
  TARGETS magic_enum
  EXPORT magic_enumTargets
  INCLUDES
  DESTINATION include
)

# Crow (uses Boost.Asio — standalone asio is not required)
CPMAddPackage(
  NAME crow
  GITHUB_REPOSITORY CrowCpp/crow
  GIT_TAG v1.2.1
  OPTIONS "CROW_BUILD_EXAMPLES OFF" "CROW_BUILD_TESTS OFF" "CROW_USE_BOOST ON"
)
target_include_directories(
  ${PROJECT_NAME} SYSTEM PUBLIC $<BUILD_INTERFACE:${crow_SOURCE_DIR}/include>
                                $<INSTALL_INTERFACE:include>
)
target_compile_definitions(${PROJECT_NAME} PUBLIC CROW_USE_BOOST)

# Torch
if(LUPNT_USE_TORCH)
  message(STATUS "[Torch] Enabled")
  find_package(Torch REQUIRED)
  target_link_libraries(${PROJECT_NAME} PUBLIC ${TORCH_LIBRARIES})
  target_include_directories(
    lupnt PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/_deps/libtorch-src/include>
                 $<INSTALL_INTERFACE:include>
  )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
else()
  message(STATUS "[Torch] Disabled")
endif()

set_target_properties(cspice PROPERTIES FOLDER thirdparty)
set_target_properties(matplot PROPERTIES FOLDER thirdparty)
set_target_properties(nodesoup PROPERTIES FOLDER thirdparty)
set_target_properties(uninstall PROPERTIES FOLDER thirdparty)

# # lupnt dependencies used for installing set(LUPNT_DEPENDENCIES "Eigen3 3.4;autodiff
# 1.1;cspice;Matplot++ 1.2;OpenMP;HighFive 2.10;Boost;yaml-cpp 0.8.0;magic_enum 0.9.7;EigenRand
# 0.5.1;OpenCV;spdlog 1.15.1" PARENT_SCOPE )

# Include PackageProject after all CPM packages are added
include(${PackageProject.cmake_SOURCE_DIR}/CMakeLists.txt)
