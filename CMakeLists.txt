# CMake project file for hibridon

cmake_minimum_required (VERSION 3.13)
project (hibridon VERSION 5.0)
enable_language (Fortran)
message("Fortran compiler ID: ${CMake_Fortran_COMPILER_ID}")
# Default build type as Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif()

# Code coverage
set(ENABLE_CODE_COVERAGE "OFF" CACHE BOOL "to enable or disable code coverage (mostly for developers)")

# profiling is disabled by default because:
# 1. most users won't need it: this is targetted to developers, not users
# 2. it makes code slower
set(ENABLE_PROFILING "OFF" CACHE BOOL "to enable or disable performance profiling (mostly for developers)")

if(ENABLE_PROFILING)
  find_program(PROFILING_GPROF_EXE gprof)
  if (NOT PROFILING_GPROF_EXE)
      message(WARNING "failed to find gprof executable (tool used to generate dpf output of performance profiling); without it, no pdf output for profiling will be produced.")
  endif()
  message(STATUS "PROFILING_GPROF_EXE = ${PROFILING_GPROF_EXE}")

  find_program(PROFILING_GPROF2DOT_EXE gprof2dot)
  if (NOT PROFILING_GPROF2DOT_EXE)
      message(WARNING "failed to find gprof2dot executable (tool used to generate dpf output of performance profiling); without it, no pdf output for profiling will be produced. gprof2dot can be installed using pip install gprof2dot.")
  endif()
  message(STATUS "PROFILING_GPROF2DOT_EXE = ${PROFILING_GPROF2DOT_EXE}")

  find_program(PROFILING_DOT_EXE dot)
  if (NOT PROFILING_DOT_EXE)
      message(WARNING "failed to find dot executable (graphviz tool used to generate dpf output of performance profiling); without it, no pdf output for profiling will be produced. dot executable is part of an open source software named graphviz.")
  endif()
  message(STATUS "PROFILING_DOT_EXE = ${PROFILING_DOT_EXE}")

  if ( PROFILING_GPROF_EXE AND PROFILING_GPROF2DOT_EXE AND PROFILING_DOT_EXE )
    set(GENERATE_PROFILING_PDF "ON")
  else()
    set(GENERATE_PROFILING_PDF "OFF")
  endif()
  message(STATUS "GENERATE_PROFILING_PDF=${GENERATE_PROFILING_PDF}")
endif()


find_package( BLAS REQUIRED )
find_package( LAPACK REQUIRED )
find_package ( OpenMP )

# get the name of the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message("Fortran compiler: ${Fortran_COMPILER_NAME}")

# Include useful functions
include(cmake/modules/add_hibexe.cmake)
include(cmake/modules/add_hibridon_test.cmake)


# Add the hibridon library
add_subdirectory(lib)

# Add the tests only if BUILD_TESTING is explicitly ON 
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
  add_subdirectory(src/unittests)
  endif()

