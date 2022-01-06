# CMake project file for hibridon

# cmake version >= 3.3 is required for IN_LIST
cmake_minimum_required (VERSION 3.3)
project (hibridon VERSION 5.0)
enable_language (Fortran)

find_package( BLAS REQUIRED )
find_package( LAPACK REQUIRED )

# Include useful functions
include(cmake/modules/add_hibexe.cmake)
include(cmake/modules/add_hibridon_test.cmake)


# Add the hibridon library
add_subdirectory(lib)

# Add the tests only if BUILD_TESTING is explicitly ON 
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()

# Add the user-supplied potentials
file (GLOB userscripts ${CMAKE_CURRENT_BINARY_DIR}/*.user.cmake)
foreach(userscript IN LISTS userscripts)
   include (${userscript})
   add_hibexe("${EXE_NAME}" "${POT_SRC_FILE}" "${p_T_MATRIX_SIZE}")
   message ("Added ${userscript} user executable")
   unset(EXE_NAME)
   unset(POT_SRC_FILE)
   unset(p_T_MATRIX_SIZE)
endforeach()