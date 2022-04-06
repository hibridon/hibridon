# CMake project file for hibridon

cmake_minimum_required (VERSION 3.13)
project (hibridon VERSION 5.0)
enable_language (Fortran)
message("Fortran compiler ID: ${CMake_Fortran_COMPILER_ID}")
# Default build type as Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package( BLAS REQUIRED )
find_package( LAPACK REQUIRED )

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

