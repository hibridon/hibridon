# CMake project file for hibridon

cmake_minimum_required (VERSION 2.6)
project (hibridon)
enable_language (Fortran)

# set(MKLROOT )
# set(BLA_VENDOR Intel10_64lp)
find_package( BLAS REQUIRED )
message( STATUS BLAS found: ${BLAS_LIBRARIES} )
find_package( LAPACK REQUIRED )
message( STATUS BLAS found: ${LAPACK_LIBRARIES} )

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

add_library(hib
     src/hiba13p.f
     src/hiba1d3p.f
     src/hiba1del.f
     src/hiba1sg1sg.f
     src/hiba1sg.f
     src/hiba22p.f
     src/hiba2del.f
     src/hiba2mol.f
     src/hiba2pi1sg.f
     src/hiba2pi.f
     src/hiba2sg1sg.f
     src/hiba2sg.f
     src/hiba3p2s.f
     src/hiba3sg1sg.f
     src/hibaastp1.f
     src/hibaastp2.f
     src/hibaastp3.f
     src/hibaastp.f
     src/hibach2x.f
     src/hibadiat2p.f
     src/hibah2p.f
     src/hibah3p.f
     src/hiamp.f
     src/hibapi.f
     src/hibasgpi1.f
     src/hibasgpi.f
     src/hibasphtp.f
     src/hibastp1.f
     src/hibastp1sg.f
     src/hibastp.f
     src/hibastpln.f
     src/hibound.f
     src/hibrid1.f
     src/hibrid2.f
     src/hibrid3.f
     src/hibrid4.f
     src/hibrid5.f
     src/hicommon.f
     src/hidifcrs.f
     src/hihypxsc.f
     src/hiiolib_f.F
     src/himatrix.f
     src/hinput.f
     src/hiprsbr.f
     src/histmix.f
     src/hisystem.f
     src/hitensor.f
     src/hitrnprt.f
     src/hiutil.F
     src/hivector.f
     ${CMAKE_CURRENT_BINARY_DIR}/hiversion.F)
target_compile_options(hib PRIVATE "-ffixed-line-length-132")

add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/SYSCONFIG
    COMMAND "PATH=${CMAKE_CURRENT_SOURCE_DIR}/bin:$ENV{PATH}" ${CMAKE_CURRENT_SOURCE_DIR}/bin/hibdetails
    COMMENT "building ${CMAKE_CURRENT_BINARY_DIR}/SYSCONFIG"
)

add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/src/common
    COMMAND mkdir -p "${CMAKE_CURRENT_BINARY_DIR}/src/common"
    COMMENT "making directory ${CMAKE_CURRENT_BINARY_DIR}/src/common"
)


add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/src/common/parhlp
    COMMAND sed "s+xxxx+${CMAKE_CURRENT_SOURCE_DIR}+" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/parhlp.t" > "${CMAKE_CURRENT_BINARY_DIR}/src/common/parhlp"
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/src/common
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/common/parhlp.t
    COMMENT "building ${CMAKE_CURRENT_BINARY_DIR}/src/common/parhlp"
)

add_custom_target(parhlp
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/src/common/parhlp
)

add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/hiversion.F
    COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/src/hiversion.t ${CMAKE_CURRENT_BINARY_DIR}/SYSCONFIG > ${CMAKE_CURRENT_BINARY_DIR}/hiversion.F
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/hiversion.t
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/SYSCONFIG
    COMMENT "building ${CMAKE_CURRENT_BINARY_DIR}/hiversion.F from ${CMAKE_CURRENT_SOURCE_DIR}/src/hiversion.t"
)

# -std=legacy so that the compiler doesn't complain when the size of arrays don't match between a function declaration and a call to this function. At the moment, there's plenty of such cases in the code; these are not bugs, just an old way of doing things.
target_compile_options(hib PRIVATE "-std=legacy")
target_compile_definitions(hib PUBLIC __LINUX__=1)
target_include_directories(hib PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
add_dependencies(hib parhlp)

add_executable(hib_arn2 src/pot/pot_arn2.f src/himain.F)
target_compile_definitions(hib_arn2 PUBLIC KMAX=151)
target_link_libraries(hib_arn2 hib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(hib_arn2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(hib_arn2 PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
