# CMake project file for hibridon

# cmake version >= 3.3 is required for IN_LIST
cmake_minimum_required (VERSION 3.3)
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
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -fsanitize=address")
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
     src/hipotutil.F
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

set(HIBRIDON_COMPILE_OPTIONS  "")

# allow for source code lines of up to 132 columns whichi is the maximum allowed in intel fortran compiler
# (true, we are compiling with gfortran here, having the same restriction here will ease detecting lines that exceed 132 columns,
# they would break intel fortran's build)
list(APPEND HIBRIDON_COMPILE_OPTIONS "-ffixed-line-length-132")

# -std=legacy so that the compiler doesn't complain when the size of arrays don't match between a function declaration and a call to this function. At the moment, there's plenty of such cases in the code; these are not bugs, just an old way of doing things.
list(APPEND HIBRIDON_COMPILE_OPTIONS "-std=legacy")

# force arrays to be allocated on the heap instead of the stack, this removes segmentation faults crash, not sure why
list(APPEND HIBRIDON_COMPILE_OPTIONS "-fmax-stack-var-size=1")

# use large model to prevent linking errors (relocation truncated to fit: R_X86_64_PC32) because the common blocks are too big to fit in the default memory model
if(UNIX AND NOT APPLE)
list(APPEND HIBRIDON_COMPILE_OPTIONS "-mcmodel=large")
endif()
#  -fbounds-check'
target_compile_options(hib PRIVATE ${HIBRIDON_COMPILE_OPTIONS})

target_compile_definitions(hib PUBLIC HIB_UNIX_X86=1)
target_compile_definitions(hib PUBLIC HIB_UNIX=1)
target_compile_definitions(hib PUBLIC HIB_UNIX_GFORTRAN=1)
target_include_directories(hib PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
add_dependencies(hib parhlp)

# list of various architectures encountered when converting cstart/cend with fpp directives
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


enable_testing()

function(add_hibridon_test TEST_ID TEST_POT_SRC_FILE TEST_POT_DATA_FILES TEST_COMMAND_FILE TEST_INPUT_FILES TEST_OUTPUT_FILES TEST_KMAX TEST_T_MATRIX_SIZE)

  set(TEST_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests/${TEST_ID})
  set(TEST_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/tests/${TEST_ID})
  set(TEST_EXE ${TEST_ID})

  add_executable(${TEST_EXE} src/pot/${TEST_POT_SRC_FILE} src/himain.F)
  target_compile_definitions(${TEST_EXE} PUBLIC KMAX=${TEST_KMAX})
  target_compile_definitions(${TEST_EXE} PUBLIC T_MATRIX_SIZE=${TEST_T_MATRIX_SIZE})
  target_link_libraries(${TEST_EXE} hib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  target_include_directories(${TEST_EXE} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_include_directories(${TEST_EXE} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
  target_compile_options(${TEST_EXE} PRIVATE ${HIBRIDON_COMPILE_OPTIONS})


  add_test(NAME ${TEST_ID}_test_setup
      COMMAND bash -c "mkdir -p ${TEST_BUILD_DIR}/potdata && \
                      for input_file in ${TEST_INPUT_FILES} ${TEST_COMMAND_FILE}; \
                      do \
                          cp ${TEST_SRC_DIR}/$input_file ${TEST_BUILD_DIR}/ ; \
                      done ; \
                      for pot_data_file in ${TEST_POT_DATA_FILES}; \
                      do \
                          cp ${TEST_SRC_DIR}/$pot_data_file ${TEST_BUILD_DIR}/potdata ; \
                      done"
  )

  add_test(NAME ${TEST_ID}_build
      COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${TEST_EXE})

  # $basedir/bin/progs/arn2_151 <arn2_test.com
  # todo: fix memory leak and remove export ASAN_OPTIONS=detect_leaks=0
  add_test(NAME ${TEST_ID}_test_run
      COMMAND bash -c "export ASAN_OPTIONS=detect_leaks=0; cat ./${TEST_COMMAND_FILE} | ${CMAKE_CURRENT_BINARY_DIR}/${TEST_EXE}"
      WORKING_DIRECTORY "${TEST_BUILD_DIR}"
  )
  if(${TEST_ID} STREQUAL "ch3i")
    # Check that outputs from test runs are the same as the ref ones
    # todo: replace with a better comparer that uses a tolerance (at the moment, we compate the meaningful text but not the last line because the numbers are so tiny that they differ from a build and another)
    add_test(NAME ${TEST_ID}_test_check
        COMMAND bash -c "mkdir -p ${TEST_BUILD_DIR}/flx_compare \
         && cat ${TEST_SRC_DIR}/Ch3itest.flx | tail -35 | head -34 > ${TEST_BUILD_DIR}/flx_compare/gold.txt \
         && cat ${TEST_BUILD_DIR}/Ch3itest.flx | tail -35 | head -34 > ${TEST_BUILD_DIR}/flx_compare/measured.txt \
         && diff ${TEST_BUILD_DIR}/flx_compare/gold.txt ${TEST_BUILD_DIR}/flx_compare/measured.txt"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
  else()

    # Check that outputs from test runs are the same as the ref ones
    add_test(NAME ${TEST_ID}_test_check
        COMMAND bash -c "set -o errexit; \
                         for output_file in ${TEST_OUTPUT_FILES} ; \
                         do \
                            ./check_outputs ${TEST_SRC_DIR}/$output_file ${TEST_BUILD_DIR}/$output_file ; 
                         done"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
  endif()

  add_test(NAME ${TEST_ID}_test_cleanup
      COMMAND rm -R ${TEST_BUILD_DIR};
  )

  set_tests_properties(${TEST_ID}_test_setup PROPERTIES FIXTURES_SETUP ${TEST_ID}_resources)
  set_tests_properties(${TEST_ID}_test_cleanup PROPERTIES FIXTURES_CLEANUP ${TEST_ID}_resources)

  set_tests_properties(${TEST_ID}_test_run PROPERTIES FIXTURES_REQUIRED ${TEST_ID}_resources)
  set_tests_properties(${TEST_ID}_test_run PROPERTIES DEPENDS ${TEST_ID}_build)
  set_tests_properties(${TEST_ID}_test_check PROPERTIES DEPENDS check_outputs_build)

  set_tests_properties(${TEST_ID}_build ${TEST_ID}_test_run ${TEST_ID}_test_check PROPERTIES LABELS ${TEST_ID})
endfunction()

# This executable compares output files from test runs and reference outputs located in the tests directory
add_executable(check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/ci/comp_tests.f90)

add_test(NAME check_outputs_build
    COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target check_outputs)


# function to declare a test suite (group of tests)
# once a test suite has beed added, the user can lanch it using : make <TESTSUITE_ID>
function(add_hibridon_testsuite TESTSUITE_ID TESTSUITE_TESTS TESTS_LIST)
  # TESTSUITE_ID : eg 'short'
  # TESTSUITE_TESTS : list of tests to include in this test suite. eg 'arn2;arno'
  # TESTS_LIST : list of declared tests (used to check if tests in testsuite have been declared)
  set(tests_regex "")
  foreach (test_id IN LISTS TESTSUITE_TESTS)
    if("${test_id}" IN_LIST TESTS_LIST)
      list(APPEND tests_regex "^${test_id}$")
    else()
      message( FATAL_ERROR "missing test ${test_id} in test suite ${TESTSUITE_ID}")
    endif()
  endforeach()

  list(JOIN tests_regex "|" testsuite_regex)

  #message(${short_testsuite})
  # message("${testsuite_regex}")
  # set(short_testsuite_regex "(arn2_test|arno_test)")
  add_custom_target(testsuite_${TESTSUITE_ID} COMMAND bash -c "${CMAKE_CTEST_COMMAND} -L '(${testsuite_regex})'" VERBATIM)
  add_dependencies(testsuite_${TESTSUITE_ID} check_outputs hib) # note: adding hib as a dependency is not mandatory, but it makes the build process of hib library more visible (warnings are visible)
endfunction()

# declare the tests

set(tests_list 
  arn2
  arn2_big
  arno
  aroh
  boh2
  ch3he_ccsdt
  ch3i
  h2ohe_orig
  heco
  hecn 
  vfit)

foreach(test_id IN LISTS tests_list)
  include("tests/${test_id}/${test_id}.cmake")
endforeach()

# declare the test suites

set(original_testsuite
  arn2
  arn2_big
  arno
  aroh
  boh2
  ch3he_ccsdt
  ch3i
  h2ohe_orig
  heco
  hecn 
  vfit)

add_hibridon_testsuite("original" "${original_testsuite}" "${tests_list}")

# test suite from autotest/testgroup_short

set(short_testsuite 
  arn2
  arno)

add_hibridon_testsuite("short" "${short_testsuite}" "${tests_list}")

# test suite from autotest/testgroup_regular

set(regular_testsuite
  arn2
  arno
#  boh2_bound
#  ch3he
#  cph2
#  h2ohe
#  hco_trnprt
  hecn
#  nh3h2
#  ohh2_j1
  vfit)

add_hibridon_testsuite("regular" "${regular_testsuite}" "${tests_list}")

# test suite from autotest/testgroup_extended

set(extended_testsuite
  arn2
  aroh
  arno
#  arno_bound
#  arch4
#  boh2_bound
#  c2hh2
#  ch3he
#  ch3he_vib
  ch3i
#  cph2
#  h2ohe
#  hco_trnprt
  hecn
  heco
#  heohx
#  heohx_ohh2pot
#  nh3h2
#  ohh2_j1
#  ohh2_j13
#  prsbr
#  stmix
#  stp1sg_pnh3ph2
#  stp1sg_och3he
#  trnprt
  vfit)

add_hibridon_testsuite("extended" "${extended_testsuite}" "${tests_list}")

# test suite from autotest/testgroup_full_qma

set(full_qma_testsuite
  arn2
  arno
#  arno_bound
  aroh
#  boh2_bound
#  ch3he
#  ch3he_vib
  ch3i
#  h2ohe
  hecn
#  heohx
#  heohx_ohh2pot
  heco
#  neoh_sumi
#  nh3h2
#  nh3h2_qma
#  ohh2_j1
#  ohh2_j13
#  ohh2_j1_bausr
#  prsbr
#  stmix
#  trnprt
  vfit)

add_hibridon_testsuite("full_qma" "${full_qma_testsuite}" "${tests_list}")

# test suite from autotest/testgroup_all

set(all_testsuite
  arn2
  aroh
  arno
#  arno_bound
  aroh
#  boh2_bound
#  ch3he
#  ch3he_vib
  ch3i
#  cph2
#  h2ohe
  hecn
  heco
#  heohx
#  heohx_ohh2pot
#  nh3h2
#  ohh2_j1
#  ohh2_j13
#  prsbr
#  stmix
#  stp1sg_pnh3ph2
#  stp1sg_och3he
#  trnprt
  vfit)

add_hibridon_testsuite("all" "${all_testsuite}" "${tests_list}")
