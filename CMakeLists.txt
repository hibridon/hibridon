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

add_library(hib
     src/assert.f90
     src/hiamp.F90
     src/hibound.F90
     src/hibrid1.F90
     src/hibrid2.F90
     src/hibrid3.F90
     src/hibrid4.F90
     src/hibrid5.F90
     src/hicommon.F90
     src/hidifcrs.F90
     src/hiflow.F90
     src/hihypxsc.F90
     src/hiiolib_f.F90
     src/himatrix.F90
     src/himem.F90
     src/hinput.F90
     src/hiprsbr.F90
     src/histmix.F90
     src/hisystem.F90
     src/hitensor.F90
     src/hitrnprt.F90
     src/hiutil.F90
     src/hivector.F90
     src/hipotutil.F90
     src/hiversion.F90
     src/bases/hiba01_1sg.F90
     src/bases/hiba02_2sg.F90
     src/bases/hiba03_2pi.F90
     src/bases/hiba04_sgpi.F90
     src/bases/hiba05_pi.F90
     src/bases/hiba06_stp.F90
     src/bases/hiba07_13p.F90
     src/bases/hiba08_2mol.F90
     src/bases/hiba09_stpln.F90
     src/bases/hiba10_22p.F90
     src/bases/hiba11_1del.F90
     src/bases/hiba12_h2p.F90
     src/bases/hiba13_h3p.F90
     src/bases/hiba14_2del.F90
     src/bases/hiba15_diat2p.F90
     src/bases/hiba16_astp.F90
     src/bases/hiba17_ch2x.F90
     src/bases/hiba18_stp1.F90
     src/bases/hiba19_sgpi1.F90
     src/bases/hiba20_2pi1sg.F90
     src/bases/hiba21_stp1sg.F90
     src/bases/hiba22_1d3p.F90
     src/bases/hiba23_3p2s.F90
     src/bases/hiba24_sphtp.F90
     src/bases/hiba25_1sg1sg.F90
     src/bases/hiba26_2sg1sg.F90
     src/bases/hiba27_astp1.F90
     src/bases/hiba28_3sg1sg.F90
     src/bases/hiba29_astp2.F90
     src/bases/hiba30_astp3.F90
     )
# So the compiler can find the src/common included files
include_directories(src)
  
# Compile options
set (HIBRIDON_COMPILE_OPTIONS "")

if (CMake_Fortran_COMPILER_ID STREQUAL "GNU")
  include(cmake/gnu_compile_options.cmake)
elseif (CMake_Fortran_COMPILER_ID STREQUAL "Intel")
  include(cmake/intel_compile_options.cmake)
endif()

target_compile_options(hib PUBLIC ${HIBRIDON_COMPILE_OPTIONS})


# The following is to define the machine type and compiler
target_compile_definitions(hib PUBLIC HIB_UNIX_X86=1)
target_compile_definitions(hib PUBLIC HIB_UNIX=1)
if (CMake_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_definitions(hib PUBLIC HIB_UNIX_GFORTRAN=1)
elseif (CMake_Fortran_COMPILER_ID STREQUAL "Intel")
  target_compile_definitions(hib PUBLIC HIB_UNIX_IFORT=1)
endif()

# The following is for preprocessing of hiversion.F90
target_compile_definitions(hib PUBLIC _BUILD_VERS_=${PROJECT_VERSION})
target_compile_definitions(hib PUBLIC _BUILD_TYPE_="${CMAKE_BUILD_TYPE}")
target_compile_definitions(hib PUBLIC _GIT_REVISION_="${GIT_REVISION}")
target_compile_definitions(hib PUBLIC _BUILD_SYS1_="${CMAKE_HOST_SYSTEM_NAME}")
target_compile_definitions(hib PUBLIC _BUILD_SYS2_="${CMAKE_HOST_SYSTEM_VERSION}")
target_compile_definitions(hib PUBLIC _BUILD_PROC_="${CMAKE_HOST_SYSTEM_PROCESSOR}")
target_compile_definitions(hib PUBLIC _BUILD_COMP_="${CMAKE_Fortran_COMPILER}")

# The ffollowing is for preprocessing of hhelp.F90 and hiutil.F90
target_compile_definitions(hib PUBLIC _HELPDIR_="${CMAKE_CURRENT_SOURCE_DIR}/doc/")
target_compile_definitions(hib PUBLIC _BUILD_FLAGS_="${hib_compile_options}")


#target_compile_definitions(hib PUBLIC _BUILD_FLAGS_="hib_compile_options")
string(TIMESTAMP TODAY "%Y-%m-%d at %H:%M")
target_compile_definitions(hib PUBLIC _BUILD_DATE_="${TODAY}")
# The ffollowing is for preprocessing of hhelp.F90 and hiutil.F90
target_compile_definitions(hib PUBLIC _HELPDIR_="${CMAKE_CURRENT_SOURCE_DIR}/doc/")
# The following is to limit the max jtot number
target_compile_definitions(hib PUBLIC MAX_NJTOT=1000)
target_compile_definitions(hib PUBLIC DISABLE_HIB_ASSERT=1)

target_link_libraries(hib INTERFACE "${LAPACK_LIBRARIES}")
target_link_libraries(hib INTERFACE "${BLAS_LIBRARIES}")

enable_testing()

# Include all cmake files from tests directories
file (GLOB_RECURSE tests RELATIVE ${hibridon_SOURCE_DIR}/tests/ *.cmake)
foreach(test IN LISTS tests)
  include(${hibridon_SOURCE_DIR}/tests/${test})
endforeach()




