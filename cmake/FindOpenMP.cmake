# This file adds OpenMP to your project if you are using Apple Clang.

find_package(OpenMP QUIET)
if(NOT ${OpenMP_FOUND} OR NOT ${OpenMP_CXX_FOUND})
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL AppleClang AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7)
    find_program(BREW NAMES brew)
    if(BREW)
      execute_process(
        COMMAND ${BREW} ls libomp
        RESULT_VARIABLE BREW_RESULT_CODE
        OUTPUT_QUIET ERROR_QUIET
      )
      if(BREW_RESULT_CODE)
        message(
          STATUS
            "This program supports OpenMP on Mac through Brew. Please run \"brew install libomp\""
        )
      else()
        execute_process(
          COMMAND ${BREW} --prefix libomp
          OUTPUT_VARIABLE BREW_LIBOMP_PREFIX
          OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        set(OpenMP_C_LIB_NAMES "libomp")
        set(OpenMP_CXX_LIB_NAMES "libomp")
        set(OpenMP_libomp_LIBRARY ${BREW_LIBOMP_PREFIX}/lib/libomp.dylib)
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp ${BREW_LIBOMP_PREFIX}/lib/libomp.dylib")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp ${BREW_LIBOMP_PREFIX}/lib/libomp.dylib")
        set(OpenMP_INCLUDE_DIRS ${BREW_LIBOMP_PREFIX}/include)
        message(STATUS "Using Homebrew libomp from ${BREW_LIBOMP_PREFIX}")
      endif()
    else()
      message(
        STATUS
          "This program supports OpenMP on Mac through Homebrew, installing Homebrew recommmended https://brew.sh"
      )
    endif()
  endif()
endif()

find_package(OpenMP REQUIRED)
if(NOT TARGET OpenMP::OpenMP_CXX)
  add_library(OpenMP_TARGET INTERFACE)
  add_library(OpenMP::OpenMP_CXX ALIAS OpenMP_TARGET)
  target_compile_options(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
  find_package(Threads REQUIRED)
  target_link_libraries(OpenMP_TARGET INTERFACE Threads::Threads)
  target_link_libraries(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
  target_include_directories(OpenMP_TARGET INTERFACE ${OpenMP_INCLUDE_DIRS})
endif()
