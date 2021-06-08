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
     src/hiba13p.F
     src/hiba1d3p.F
     src/hiba1del.F
     src/hiba1sg1sg.F
     src/hiba1sg.F
     src/hiba22p.F
     src/hiba2del.F
     src/hiba2mol.F
     src/hiba2pi1sg.F
     src/hiba2pi.F
     src/hiba2sg1sg.F
     src/hiba2sg.F
     src/hiba3p2s.F
     src/hiba3sg1sg.F
     src/hibaastp1.F
     src/hibaastp2.F
     src/hibaastp3.F
     src/hibaastp.F
     src/hibach2x.F
     src/hibadiat2p.F
     src/hibah2p.F
     src/hibah3p.F
     src/hiamp.F
     src/hibapi.F
     src/hibasgpi1.F
     src/hibasgpi.F
     src/hibasphtp.F
     src/hibastp1.F
     src/hibastp1sg.F
     src/hibastp.F
     src/hibastpln.F
     src/hibound.F
     src/hibrid1.F
     src/hibrid2.F
     src/hibrid3.F
     src/hibrid4.F
     src/hibrid5.F
     src/hicommon.F
     src/hidifcrs.F
     src/hihypxsc.F
     src/hiiolib_f.F
     src/himatrix.F
     src/hinput.F
     src/hiprsbr.F
     src/histmix.F
     src/hisystem.F
     src/hitensor.F
     src/hitrnprt.F
     src/hiutil.F
     src/hivector.F
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

# allow for source code lines of up to 132 columns whichi is the maximum allowed in intel fortran compiler
# (true, we are compiling with gfortran here, having the same restriction here will ease detecting lines that exceed 132 columns,
# they would break intel fortran's build)
target_compile_options(hib PRIVATE "-ffixed-line-length-132")

# -std=legacy so that the compiler doesn't complain when the size of arrays don't match between a function declaration and a call to this function. At the moment, there's plenty of such cases in the code; these are not bugs, just an old way of doing things.
target_compile_options(hib PRIVATE "-std=legacy")

# force arrays to be allocated on the heap instead of the stack, this removes segmentation faults crash, not sure why
target_compile_options(hib PRIVATE "-fmax-stack-var-size=1")

# use large model to prevent linking errors (relocation truncated to fit: R_X86_64_PC32) because the common blocks are too big to fit in the default memory model
target_compile_options(hib PRIVATE "-mcmodel=large")

#  -fbounds-check'

target_compile_definitions(hib PUBLIC HIB_UNIX_X86=1)
target_compile_definitions(hib PUBLIC HIB_UNIX=1)
target_compile_definitions(hib PUBLIC HIB_UNIX_GFORTRAN=1)
target_include_directories(hib PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
add_dependencies(hib parhlp)

# BATCH
# CRAY
# CRAY_COS
# CRAY_UNICOS
# DEC_RISC
# DISCO_F
# DISCO_U
# ETA
# FPS
# IBM
# IBM_RISC
# IBM_VM
# MOLPRO
# NONE
# ULTRIX_DEC
# UNIVAC
# UNIX
# UNIX_AIX
# UNIX_BLAS3
# UNIX_CONVEX
# UNIX_DEC
# UNIX_GFORTRAN
# UNIX_HP
# UNIX_IBM
# UNIX_IFORT
# UNIX_IRIS
# UNIX_IRIX
# UNIX_MAC
# UNIX_NEC
# UNIX_NOBLAS
# UNIX_PGI
# UNIX_S
# UNIX_SEQUENT
# UNIX_SUN
# UNIX_X86
# UNIX_XLF
# UNKNOWN
# VAX

add_executable(hib_arn2 src/pot/pot_arn2.F src/himain.F)
target_compile_definitions(hib_arn2 PUBLIC KMAX=151)
target_link_libraries(hib_arn2 hib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(hib_arn2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(hib_arn2 PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)

# allow for source code lines of up to 132 columns whichi is the maximum allowed in intel fortran compiler
# (true, we are compiling with gfortran here, having the same restriction here will ease detecting lines that exceed 132 columns,
# they would break intel fortran's build)
target_compile_options(hib_arn2 PRIVATE "-ffixed-line-length-132")

# -std=legacy so that the compiler doesn't complain when the size of arrays don't match between a function declaration and a call to this function. At the moment, there's plenty of such cases in the code; these are not bugs, just an old way of doing things.
target_compile_options(hib_arn2 PRIVATE "-std=legacy")

# force arrays to be allocated on the heap instead of the stack, this removes segmentation faults crash, not sure why
target_compile_options(hib_arn2 PRIVATE "-fmax-stack-var-size=1")

# use large model to prevent linking errors (relocation truncated to fit: R_X86_64_PC32) because the common blocks are too big to fit in the default memory model
target_compile_options(hib_arn2 PRIVATE "-mcmodel=large")


enable_testing()

# This executable compares output files from test runs and reference outputs located in the tests directory
add_executable(check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/ci/comp_tests.f90)

add_test(NAME check_outputs_build
    COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target check_outputs)

add_test(NAME arn2_test_setup
    COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/arn2; \
                    cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/arn2_test.com \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Arn2_test.inp \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Arn2_dxsec.inp \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Arn2.fluxinp \
                        ${CMAKE_CURRENT_BINARY_DIR}/tests/arn2/"
)

add_test(NAME arn2_build
    COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hib_arn2)

# $basedir/bin/progs/hib_arn2_151 <arn2_test.com
add_test(NAME arn2_test_run
    COMMAND bash -c "cat ./arn2_test.com | ${CMAKE_CURRENT_BINARY_DIR}/hib_arn2"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests/arn2"
)

# Check that outputs from test runs are the same as the ref ones
add_test(NAME arn2_test_check
    COMMAND bash -c "./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Cctest1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/arn2/Cctest1.ics && \
    ./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Ccrstest1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/arn2/Ccrstest1.ics && \
    ./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Cstest1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/arn2/Cstest1.ics"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_test(NAME arn2_test_cleanup
    COMMAND ls -R ${CMAKE_CURRENT_BINARY_DIR}/tests/arn2;
)


set_tests_properties(arn2_test_setup PROPERTIES FIXTURES_SETUP arn2)
set_tests_properties(arn2_test_cleanup PROPERTIES FIXTURES_CLEANUP arn2)

set_tests_properties(arn2_test_run PROPERTIES FIXTURES_REQUIRED arn2)
set_tests_properties(arn2_test_run PROPERTIES DEPENDS arn2_build)
set_tests_properties(arn2_test_check PROPERTIES DEPENDS check_outputs_build)
