# CMake project file for hibridon

# cmake version >= 3.3 is required for IN_LIST
cmake_minimum_required (VERSION 3.3)
project (hibridon VERSION 5.0)
enable_language (Fortran)

# code coverage is disabled by default because:
# 1. most users won't need it: this is targetted to developers, not users
# 2. code coverage requires extra programs, and this would make hibridon's build process more complicate than needed for normal users
set(ENABLE_CODE_COVERAGE "OFF" CACHE BOOL "to enable or disable code coverage analysis (mostly for developers)")

# profiling is disabled by default because:
# 1. most users won't need it: this is targetted to developers, not users
# 2. it makes code slower
set(ENABLE_PROFILING "OFF" CACHE BOOL "to enable or disable performance profiling (mostly for developers)")

set(ENABLE_ASSERTS "ON" CACHE BOOL "to enable or disable assert checks for debug builds (mostly for developers).")

# With gfortran, variables are initialized to 0 by default, but it's not the case with ifort compiler. Therefore, all variables need to be initialized to guarantee the program works in a deterministic way on every compiler.
set(ENABLE_UNINIT_VAR_RUNTIME_DETECTOR "ON" CACHE BOOL "if set, causes the use of uninitilazed variables to trigger an exception, thus making the program stop. This is useful to detect uninitialized memory issues.")

# store variables that are used in functions (add_hibexe and add_hibridon_test) in the cache so that external cmake projects can use these functions (otherwise these variables would be empty when used from other cmake projects)

set(HIBRIDON_COMPILE_OPTIONS "" CACHE STRING "list of compilation options for hibridon source code") 

set(HIBRIDON_LINK_LIBRARIES "" CACHE STRING "list of libraries that need to be linked with hibridon") 

set(GCC_COVERAGE_LINK_FLAGS "" CACHE STRING "the libraries that need to be linked with hibridon for code coverage") 

if(ENABLE_CODE_COVERAGE)
  find_program(COVERAGE_LCOV_EXE lcov)
  if (NOT COVERAGE_LCOV_EXE)
      message(FATAL_ERROR "failed to find lcov executable (tool used incode coverage analysis)")
  endif()
  message(STATUS "COVERAGE_LCOV_EXE = ${COVERAGE_LCOV_EXE}")

  find_program(COVERAGE_GENHTML_EXE genhtml)
  if (NOT COVERAGE_GENHTML_EXE)
      message(FATAL_ERROR "failed to find genhtml executable (tool used in code coverage analysis)")
  endif()
  message(STATUS "COVERAGE_GENHTML_EXE = ${COVERAGE_GENHTML_EXE}")
endif()

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

endif()

# set(MKLROOT )
# set(BLA_VENDOR Intel10_64lp)
find_package( BLAS REQUIRED )
message( STATUS BLAS found: ${BLAS_LIBRARIES} )
find_package( LAPACK REQUIRED )
message( STATUS BLAS found: ${LAPACK_LIBRARIES} )

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
     src/bases/hiba30_astp3.F90)


# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")

  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -march=native")
  set (warning_options "-Wall -Wextra") # activate all wanings

  set (warning_options "${warning_options} -Wno-do-subscript")
  set (warning_options "${warning_options} -Wno-unused-dummy-argument")
  set (warning_options "${warning_options} -Wno-unused-variable")
  set (warning_options "${warning_options} -Wno-unused-parameter")
  set (warning_options "${warning_options} -Wno-conversion")
  set (warning_options "${warning_options} -Wno-argument-mismatch")
  set (warning_options "${warning_options} -Wno-unused-label")
  set (warning_options "${warning_options} -Wno-tabs")
  set (warning_options "${warning_options} -Wno-compare-reals")
  set (warning_options "${warning_options} -Wno-intrinsic-shadow")
  set (warning_options "${warning_options} -Wno-character-truncation")
  

  # don't initialize local variables to zero : initialize them with dummy values to force failures which help detect uninitialized variables as -Wuninitialized doesn't detect everything (see issue #38 or https://stackoverflow.com/questions/39591300/gfortran-wuninitialized-flag-misses-variable-in-do-loop)
  # by default, gfortran initializes integers to 0, but not ifort : as a result, some bugs in the code are hidden with gfortran default options
  # set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fno-init-local-zero")
  if(ENABLE_UNINIT_VAR_RUNTIME_DETECTOR)
    # initialize variables to something else than 0 to force the programe to behave badly in case of unitialized variables
    set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -finit-integer=333333333")
    set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -finit-real=snan")
    set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ffpe-trap=invalid,zero,overflow")
    target_compile_definitions(hib PUBLIC FLOATING_POINT_ERRORS_CAUSE_EXCEPTIONS=1)
  else()
    # from gfortran documentation :
    # Finally, note that enabling any of the -finit-* options will silence warnings that would have been emitted by -Wuninitialized for the affected local variables.
    set (warning_options "${warning_options} -Wuninitialized")
  endif()

  set (CMAKE_Fortran_FLAGS_DEBUG "-fno-f2c -O0 -g  -fsanitize=address ${warning_options}")

  # todo: add -fbounds-check once all mismatched declarations of variables (eg in common blocks) have been fixed
  if(ENABLE_CODE_COVERAGE)
    set(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage --coverage")
    set(GCC_COVERAGE_LINK_FLAGS "-lgcov")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${GCC_COVERAGE_COMPILE_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${GCC_COVERAGE_LINK_FLAGS}")
  endif()
  if(ENABLE_PROFILING)
    set(GCC_PROFILING_COMPILE_FLAGS "-pg")
    set(GCC_PROFILING_LINK_FLAGS "-pg")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${GCC_PROFILING_COMPILE_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${GCC_PROFILING_LINK_FLAGS}")
  endif()
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort

  set (warning_options "")
  # set (warning_options "-warn all") # activate all warnings

  # [no]alignments
  #     Determines whether warnings occur for data that is not naturally aligned. 
  set (warning_options "${warning_options} -warn alignments")

  # [no]declarations
  #     Determines whether warnings occur for any undeclared names. 
  set (warning_options "${warning_options} -warn declarations")

  # [no]errors
  #     Determines whether warnings are changed to errors. 
  # set (warning_options "${warning_options} -warn noerrors")

  # [no]externals
  #     Determines whether warnings occur for any dummy procedures or procedure calls that have no explicit interface or have not been declared external. 
  set (warning_options "${warning_options} -warn externals")

  # [no]general
  #     Determines whether warning messages and informational messages are issued by the compiler. 
  set (warning_options "${warning_options} -warn general")

  # [no]ignore_loc
  #     Determines whether warnings occur when %LOC is stripped from an actual argument. 
  set (warning_options "${warning_options} -warn ignore_loc")

  # [no]interfaces
  #     Determines whether the compiler checks the interfaces of all SUBROUTINEs called and FUNCTIONs invoked in your compilation against an external set of interface blocks.
  # this option causes the following errors (that can't be disabled by diag-disable):
  # error #6633: The type of the actual argument differs from the type of the dummy argument.
  # error #6634: The shape matching rules of actual arguments and dummy arguments have been violated.
  # error #7983: The storage extent of the dummy argument exceeds that of the actual argument.
  # set (warning_options "${warning_options} -warn interfaces")


  # [no]shape
  #     Determines whether array conformance violations are diagnosed with errors or warnings when used with the check shape
  #     option. 
  set (warning_options "${warning_options} -warn shape")

  # [no]stderrors
  #     Determines whether warnings about Fortran standard violations are changed to errors.

  # [no]truncated_source
  #     Determines whether warnings occur when source exceeds the maximum column width in fixed-format files. 
  set (warning_options "${warning_options} -warn truncated_source")

  # [no]uncalled
  #     Determines whether warnings occur when a statement function is never called 
  set (warning_options "${warning_options} -warn uncalled")

  # [no]unused
  #     Determines whether warnings occur for declared variables that are never used. 
  set (warning_options "${warning_options} -warn unused")

  # [no]usage
  #     Determines whether warnings occur for questionable programming practices. 
  set (warning_options "${warning_options} -warn usage")

  # disable warning #6717: This name has not been given an explicit type
  set (warning_options "${warning_options} -diag-disable 6717")

  # disable remark #7712: This variable has not been used.
  set (warning_options "${warning_options} -diag-disable 7712")

  # disable warning #8889: Explicit declaration of the EXTERNAL attribute is required.
  set (warning_options "${warning_options} -diag-disable 8889")
  

  set (common_options "")

  if(ENABLE_UNINIT_VAR_RUNTIME_DETECTOR)
    # Force à NaN toutes les variables de type intrinsèque ainsi que les tableaux non initialisés. Cette option implique -fpe0. Pour éviter de récupérer des exceptions, qui ne soient pas relatives à des variables non-initialisées, nous recommandons de réduire le niveau d'optimisation à -O1 ou -O0 ou alors d'utiliser -fp-speculation=safe pour faire la détection de variables non-initialisés.
    set (common_options "${common_options} -init=arrays,snan")
  endif()

  # Cette combinaison d'options stoppe l'exécution dès qu'une exception (overflow, underflow, division par zéro, opération invalide,…) se produit; elle indique à quel niveau du code elle s'est produite. Le contrôle s'opère dans chaque subroutine, contrairement à l'option -fpe0 qui agit uniquement sur le programme principal. 
#  set (common_options "${common_options} -fpe-all=0 -no-ftz -traceback")

  # set (warning_options "") # activate all warnings
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${common_options} ${warning_options}")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback -fp-stack-check ${common_options} ${warning_options}")

  if(ENABLE_UNINIT_VAR_RUNTIME_DETECTOR)
    # check uses of uninitialized variables in run time
    # this is very useful as these are always bugs that cause the program to behave randomly with ifort (gfortran initializes data with zero)
    set (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -check uninit")
    target_compile_definitions(hib PUBLIC FLOATING_POINT_ERRORS_CAUSE_EXCEPTIONS=1)
  endif()


endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")




target_compile_definitions(hib PUBLIC MAX_NJTOT=1000)

set(GIT_REVISION 'not-available')
execute_process(COMMAND "git" "--git-dir=${CMAKE_CURRENT_SOURCE_DIR}/.git" "log" "--pretty=format:%h" "-n" "1"
                OUTPUT_VARIABLE GIT_REVISION
                ERROR_QUIET)
                

# The following is for preprocessing of hiversion.F90
target_compile_definitions(hib PUBLIC _BUILD_VERS_=${PROJECT_VERSION})
target_compile_definitions(hib PUBLIC _BUILD_TYPE_="${CMAKE_BUILD_TYPE}")
target_compile_definitions(hib PUBLIC _GIT_REVISION_="${GIT_REVISION}")
target_compile_definitions(hib PUBLIC _BUILD_SYS1_="${CMAKE_HOST_SYSTEM_NAME}")
target_compile_definitions(hib PUBLIC _BUILD_SYS2_="${CMAKE_HOST_SYSTEM_VERSION}")
target_compile_definitions(hib PUBLIC _BUILD_PROC_="${CMAKE_HOST_SYSTEM_PROCESSOR}")
target_compile_definitions(hib PUBLIC _BUILD_COMP_="${CMAKE_Fortran_COMPILER}")
target_compile_definitions(hib PUBLIC _BUILD_FLAGS_="${HIBRIDON_COMPILE_OPTIONS}")
string(TIMESTAMP TODAY "%Y-%m-%d at %H:%M")
target_compile_definitions(hib PUBLIC _BUILD_DATE_="${TODAY}")
# The ffollowing is for preprocessing of hhelp.F90 and hiutil.F90
target_compile_definitions(hib PUBLIC _HELPDIR_="${CMAKE_CURRENT_SOURCE_DIR}/doc/")

add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/src/common
    COMMAND mkdir -p "${CMAKE_CURRENT_BINARY_DIR}/src/common"
    COMMENT "making directory ${CMAKE_CURRENT_BINARY_DIR}/src/common"
)


set(HIBRIDON_COMPILE_OPTIONS "") # intended to be a list, not a string

list(APPEND HIBRIDON_COMPILE_OPTIONS "-I${CMAKE_CURRENT_SOURCE_DIR}/src")

# allow for source code lines of up to 132 columns which is the maximum allowed in intel fortran compiler
# (true, we are compiling with gfortran here, having the same restriction here will ease detecting lines that exceed 132 columns,
# they would break intel fortran's build)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  list(APPEND HIBRIDON_COMPILE_OPTIONS "-ffree-line-length-none")
  list(APPEND HIBRIDON_COMPILE_OPTIONS "-ffixed-line-length-132")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  list(APPEND HIBRIDON_COMPILE_OPTIONS "-extend-source" "132")
  # source code lines can be of any length, as compilation options are included as source code strings (see bug #5)
  list(APPEND HIBRIDON_COMPILE_OPTIONS "-no-wrap-margin")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")


if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # -std=legacy so that the compiler doesn't complain when the size of arrays don't match between a function declaration and a call to this function. At the moment, there's plenty of such cases in the code; these are not bugs, just an old way of doing things.
  list(APPEND HIBRIDON_COMPILE_OPTIONS "-std=legacy")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")


if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # force arrays to be allocated on the heap instead of the stack, this removes segmentation faults crash, not sure why
  list(APPEND HIBRIDON_COMPILE_OPTIONS "-fmax-stack-var-size=1")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # force arrays to be allocated on the heap instead of the stack, this removes segmentation faults crash, not sure why
  # the following option means : automatic arrays and arrays created for temporary computations are allocated on the stack if their size can be determined at compile time and if it doesn't exceed 10kb.
  # this option seems necessary for big values of kmax (eg kmax=5000), according to bdesrousseaux, otherwise the user will experience a segmentation fault. I guess that without this option, the stack becomes too small to contain such big arrays...
  list(APPEND HIBRIDON_COMPILE_OPTIONS "-heap-arrays 10")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# use large model to prevent linking errors (relocation truncated to fit: R_X86_64_PC32) because the common blocks are too big to fit in the default memory model
if(UNIX AND NOT APPLE)
    list(APPEND HIBRIDON_COMPILE_OPTIONS "-mcmodel=large")
endif()

if(CMAKE_BUILD_TYPE MATCHES Release)
  target_compile_definitions(hib PUBLIC DISABLE_HIB_ASSERT=1)
else()
  # asserts can be activated on non-release builds (eg debug build)
  if(NOT ENABLE_ASSERTS)
    target_compile_definitions(hib PUBLIC DISABLE_HIB_ASSERT=1)
  endif()
endif()

set(HIBRIDON_LINK_LIBRARIES "")
list(APPEND HIBRIDON_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
list(APPEND HIBRIDON_LINK_LIBRARIES "${BLAS_LIBRARIES}")




#  -fbounds-check'
target_compile_options(hib PRIVATE ${HIBRIDON_COMPILE_OPTIONS})

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  target_compile_definitions(hib PUBLIC HIB_UNIX_X86=1)
  target_compile_definitions(hib PUBLIC HIB_UNIX=1)
  target_compile_definitions(hib PUBLIC HIB_UNIX_GFORTRAN=1)
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  target_compile_definitions(hib PUBLIC HIB_UNIX_X86=1)
  target_compile_definitions(hib PUBLIC HIB_UNIX=1)
  target_compile_definitions(hib PUBLIC HIB_UNIX_IFORT=1)
endif()


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

function(add_hibridon_test TEST_ID TEST_POT_SRC_FILE TEST_POT_DATA_FILES TEST_COMMAND_FILE TEST_INPUT_FILES TEST_OUTPUT_FILES TEST_T_MATRIX_SIZE)
  # message("for test ${TEST_ID}, CMAKE_CURRENT_LIST_DIR=${CMAKE_CURRENT_LIST_DIR}")
  get_filename_component(TEST_DIRNAME ${CMAKE_CURRENT_LIST_DIR} NAME)

  set(TEST_SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
  set(TEST_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/tests/${TEST_ID})
  set(TEST_EXE ${TEST_ID})
  set(HIBRIDON_MODULES_DIR ${hibridon_BINARY_DIR}) # location where to find the module files (*.mod)

  add_executable(${TEST_EXE} ${TEST_SRC_DIR}/${TEST_POT_SRC_FILE} ${hibridon_SOURCE_DIR}/src/himain.F90)
  target_compile_definitions(${TEST_EXE} PUBLIC T_MATRIX_SIZE=${TEST_T_MATRIX_SIZE})
  target_link_libraries(${TEST_EXE} hib ${HIBRIDON_LINK_LIBRARIES})
  target_include_directories(${TEST_EXE} PUBLIC ${hibridon_SOURCE_DIR}/src)
  target_include_directories(${TEST_EXE} PUBLIC ${HIBRIDON_MODULES_DIR}) 
  target_compile_options(${TEST_EXE} PRIVATE ${HIBRIDON_COMPILE_OPTIONS})

  if(ENABLE_CODE_COVERAGE)
    target_link_options(${TEST_EXE} PUBLIC ${GCC_COVERAGE_LINK_FLAGS})
  endif()


  add_test(NAME ${TEST_ID}_test_setup
      COMMAND bash -c "mkdir -p ${TEST_BUILD_DIR}/potdata && \
                      mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/coverage && \
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

  set(test_command_prefix "")

  # set(mem_leak_disabled_tests arn2 arn2_big)
  # if(${TEST_ID} IN_LIST mem_leak_disabled_tests)
  #   message( WARNING "memory leak detection has been disabled on test ${TEST_ID} until its leak is fixed")
  #   set(test_command_prefix "export ASAN_OPTIONS=detect_leaks=0; ")
  #   # todo: fix memory leak and remove export ASAN_OPTIONS=detect_leaks=0
  # endif()

  add_test(NAME ${TEST_ID}_test_run
    COMMAND bash -c "${test_command_prefix} cat ./${TEST_COMMAND_FILE} | ${CMAKE_CURRENT_BINARY_DIR}/${TEST_EXE} --kmax ${TEST_KMAX} | tee ./${TEST_ID}.stdout"
    WORKING_DIRECTORY "${TEST_BUILD_DIR}"
  )

    # Check that outputs from test runs are the same as the ref ones
  add_test(NAME ${TEST_ID}_test_check
      COMMAND bash -c "set -o errexit; \
                        for output_file in ${TEST_OUTPUT_FILES} ; \
                        do \
                          ${hibridon_BINARY_DIR}/check_outputs ${TEST_SRC_DIR}/$output_file ${TEST_BUILD_DIR}/$output_file ; 
                        done"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )

  if(ENABLE_CODE_COVERAGE)
    add_test(NAME ${TEST_ID}_save_coverage
      COMMAND lcov -c -d ${CMAKE_CURRENT_BINARY_DIR}/ -o ${CMAKE_CURRENT_BINARY_DIR}/coverage/${TEST_ID}.info;
    )

    add_test(NAME ${TEST_ID}_cleanup_coverage
      COMMAND find ${CMAKE_CURRENT_BINARY_DIR} -name "*.gcda" -delete;
    )
  endif()

  if(GENERATE_PROFILING_PDF)
    add_test(NAME ${TEST_ID}_build_profiling_pdf
      COMMAND  bash -c "${PROFILING_GPROF_EXE} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_EXE}  | ${PROFILING_GPROF2DOT_EXE} | ${PROFILING_DOT_EXE} -Tpdf -o call_graph.pdf"
      WORKING_DIRECTORY "${TEST_BUILD_DIR}"
    )
  endif()

  add_test(NAME ${TEST_ID}_test_cleanup
      COMMAND ls -R ${TEST_BUILD_DIR};
  )

  set_tests_properties(${TEST_ID}_test_setup PROPERTIES FIXTURES_SETUP ${TEST_ID}_resources)
  if(ENABLE_CODE_COVERAGE)
    set_tests_properties(${TEST_ID}_cleanup_coverage PROPERTIES FIXTURES_CLEANUP ${TEST_ID}_resources)
  endif()
  set_tests_properties(${TEST_ID}_test_cleanup PROPERTIES FIXTURES_CLEANUP ${TEST_ID}_resources)

  set_tests_properties(${TEST_ID}_test_run PROPERTIES FIXTURES_REQUIRED ${TEST_ID}_resources)
  set_tests_properties(${TEST_ID}_test_run PROPERTIES DEPENDS ${TEST_ID}_build)
  set_tests_properties(${TEST_ID}_test_check PROPERTIES DEPENDS check_outputs_build)

  if(ENABLE_CODE_COVERAGE)
    set_tests_properties(${TEST_ID}_save_coverage PROPERTIES DEPENDS ${TEST_ID}_test_check)
    set_tests_properties(${TEST_ID}_cleanup_coverage PROPERTIES DEPENDS ${TEST_ID}_save_coverage)
  endif()

  set_tests_properties(${TEST_ID}_build ${TEST_ID}_test_run ${TEST_ID}_test_check PROPERTIES LABELS ${TEST_ID})
  if(ENABLE_CODE_COVERAGE)
    set_tests_properties(${TEST_ID}_save_coverage ${TEST_ID}_cleanup_coverage PROPERTIES LABELS ${TEST_ID})
  endif()
  if(GENERATE_PROFILING_PDF)
    set_tests_properties(${TEST_ID}_build_profiling_pdf PROPERTIES LABELS ${TEST_ID})
  endif()

endfunction()

# This executable compares output files from test runs and reference outputs located in the tests directory
add_executable(check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/ci/comp_mods.f90 ${CMAKE_CURRENT_SOURCE_DIR}/ci/comp_tests.f90)

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
    if(TARGET "${test_id}")
      list(APPEND tests_regex "^${test_id}$")
    else()
      message( FATAL_ERROR "missing test ${test_id} in test suite ${TESTSUITE_ID}")
    endif()
  endforeach()

  # build testsuite_regex, which is expected to contain something like "^arn2$|^arn2_big$|^arno$|^aroh$|^boh2_b$|^ch3he_ccsdt$|^ch3i$|^h2ohe_orig$|^hecn$|^heco$|^vfit$"
  if(${CMAKE_VERSION} VERSION_LESS "3.12.0")
    set(testsuite_regex_tmp "")
    foreach(test_regex ${tests_regex})
      string(APPEND testsuite_regex_tmp "${test_regex}|")
    endforeach()
    string(REGEX REPLACE "\\|$" "" testsuite_regex "${testsuite_regex_tmp}")
  else()
    # join does exactly what we need here, but required recent versions of cmake
    list(JOIN tests_regex "|" testsuite_regex)
  endif()

  #message(${short_testsuite})
  # message("${testsuite_regex}")
  # set(short_testsuite_regex "(arn2_test|arno_test)")
  add_custom_target(testsuite_${TESTSUITE_ID} COMMAND bash -c "${CMAKE_CTEST_COMMAND} -L '(${testsuite_regex})'" VERBATIM)
  add_dependencies(testsuite_${TESTSUITE_ID} check_outputs hib) # note: adding hib as a dependency is not mandatory, but it makes the build process of hib library more visible (warnings are visible)
endfunction()

if(ENABLE_CODE_COVERAGE)
  # build 
  add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/coverage/total.info COMMAND bash -c "\
    command='${COVERAGE_LCOV_EXE} ' ; \
    for test_coverage in $(find ${CMAKE_BINARY_DIR} -name \"*.info\") ; \
    do \
      echo adding test coverage $test_coverage to ${CMAKE_BINARY_DIR}/coverage/total.info ; \
      command=\"\$command -a \$test_coverage\" ; \
    done ; \
    command=\"\$command -o ${CMAKE_BINARY_DIR}/coverage/total.info\" ; \
    \$command" VERBATIM
    COMMENT "building total code coverage ${CMAKE_BINARY_DIR}/coverage/total.info from all code coverage files found in $(find ${CMAKE_BINARY_DIR}")

  add_custom_target(html_coverages DEPENDS ${CMAKE_BINARY_DIR}/coverage/total.info VERBATIM COMMAND bash -c "\
    for test_coverage_file in $(find ${CMAKE_BINARY_DIR}/coverage -name \"*.info\") ; \
    do \
      html_root=${CMAKE_BINARY_DIR}/coverage/\$(basename \$test_coverage_file | sed 's/.info//') ; \
      echo building $html_root/index.html from $test_coverage_file ; \
      mkdir -p \$html_root ; \
      echo ${COVERAGE_GENHTML_EXE} ; \
      ${COVERAGE_GENHTML_EXE} \$test_coverage_file -o \$html_root > /dev/null ; \
    done ;"
    COMMENT "building an html report for each code coverage file found in ${CMAKE_BINARY_DIR}/coverage")
endif()

# declare the tests

set(test_dirs_list 
  arch4
  arn2
  arn2_big
  arno
  arno_bound
  aroh
  boh2
  c2hh2
  cahe
  ch2he
  ch3he
  ch3he_ccsdt
  ch3he_vib
  ch3i
  cnph2
  cph2
  h2ohe
  h2ohe_orig
  h2ohe_prsbr
  hco_trnprt 
  hecn 
  heco
  heohx
  nh3h2
  nh3h2_qma
  nh3he
  nhph2
  o3ph2
  ohh2_bausr
  ohh2_j1
  ohh2_j13
  stp1sg
  vfit)

foreach(test_dir IN LISTS test_dirs_list)
  include("tests/${test_dir}/${test_dir}.cmake")
endforeach()

# declare the test suites

set(original_testsuite
  arn2
  arn2_big
  arno
  aroh
  boh2_b
  ch3he_ccsdt
  ch3i
  h2ohe_orig
  hecn 
  heco
  vfit)

add_hibridon_testsuite("original" "${original_testsuite}" "${test_dirs_list}")

# test suite from autotest/testgroup_short

set(short_testsuite 
  arn2
  arno)

add_hibridon_testsuite("short" "${short_testsuite}" "${test_dirs_list}")

# test suite from autotest/testgroup_regular

set(regular_testsuite
  arn2
  arno
  boh2_b2
  cahe
  ch3he
  cph2
  h2ohe
  hco_trnprt
  hecn
  nh3h2
  ohh2_j1
  vfit)

add_hibridon_testsuite("regular" "${regular_testsuite}" "${test_dirs_list}")

# test suite from autotest/testgroup_extended

set(extended_testsuite
  arch4
  arn2
  arno
  arno_bound
  aroh
  boh2_b2
  c2hh2
  cahe
  ch2he_a
  ch2he_x
  ch2he_stmix
  ch3he
  ch3he_vib
  ch3i
  cph2
  h2ohe
  h2ohe_prsbr
  hco_trnprt
  hecn
  heco
  heohx
  nh3h2
  ohh2_bausr_heohx
  ohh2_j1
  ohh2_j13
  stp1sg_ch3he
  stp1sg_nh3h2
  vfit)

add_hibridon_testsuite("extended" "${extended_testsuite}" "${test_dirs_list}")

# test suite from autotest/testgroup_full_qma

set(full_qma_testsuite
  arn2
  arno
  arno_bound
  aroh
  boh2_b2
  cahe
  ch2he_a
  ch2he_x
  ch2he_stmix
  ch3he
  ch3he_vib
  ch3i
  h2ohe
  h2ohe_prsbr
  h2ohe_trnprt
  hecn
  heco
  heohx
  nh3h2
  nh3h2_qma
  ohh2_bausr_heohx
  ohh2_bausr_j1
  ohh2_j1
  ohh2_j13
  vfit)

add_hibridon_testsuite("full_qma" "${full_qma_testsuite}" "${test_dirs_list}")

# test suite from autotest/testgroup_all

set(all_testsuite
  arch4
  arch4_quick
  arn2
  arn2_big
  arno
  arno_bound
  aroh
  boh2_b
  boh2_b2
  c2hh2
  c2hh2_quick
  cahe
  ch2he_a
  ch2he_x
  ch2he_x_quick
  ch2he_stmix
  ch3he
  ch3he_quick
  ch3he_ccsdt
  ch3he_ccsdt_quick
  ch3he_vib
  ch3i
  cnph2
  cph2
  h2ohe
  h2ohe_orig
  h2ohe_prsbr
  h2ohe_prsbr_quick
  h2ohe_trnprt
  hco_trnprt
  hco_trnprt_quick
  hecn
  heco
  heohx
  nh3h2
  nh3h2_qma
  nhph2
  nhph2_quick
  o3ph2
  o3ph2_quick
  ohh2_bausr_heohx
  ohh2_bausr_j1
  ohh2_j1
  ohh2_j13
  stp1sg_ch3he
  stp1sg_nh3h2
  vfit)


add_hibridon_testsuite("all" "${all_testsuite}" "${test_dirs_list}")

# a suite which goal is to cover as much as possible hibridon's code,
# ignoring the meaningfulness of the tests
# this test suite is the one to use to detect possible regressions after refactoring code
set(coverage_testsuite
  arch4_quick
  arn2
  arn2_big
  arno
  arno_bound
  aroh
  boh2_b
  boh2_b2
  c2hh2_quick
  cahe
  ch2he_a_quick
  ch2he_x_quick
  ch2he_stmix
  ch3he_quick
  ch3he_ccsdt_quick
  ch3he_vib_quick
  ch3i
  cnph2
  h2ohe
  h2ohe_orig
  h2ohe_prsbr_quick
  h2ohe_trnprt
  hco_trnprt_quick
  hecn_quick
  heco
  heohx_quick
  nh3h2
  nh3h2_qma
  nh3he
  nhph2_quick
  o3ph2_quick
  ohh2_bausr_heohx_quick
  ohh2_bausr_j1_quick
  ohh2_j1_quick
  ohh2_j13_quick
  stp1sg_ch3he_quick
  stp1sg_nh3h2
  vfit)
# test cph2 removed until it is fixed (see issue #48)


add_hibridon_testsuite("coverage" "${coverage_testsuite}" "${test_dirs_list}")

# add_hibridon_testsuite("empty" "" "${test_dirs_list}")

# The following is to build hibridon with an user-supplied potential

function(add_hibexe EXE_NAME POT_SRC_FILE p_T_MATRIX_SIZE)
  message("Adding executable ${EXE_NAME} using potential ${POT_SRC_FILE}")
  set(TEST_SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
  set(TEST_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/tests/${TEST_ID})
  set(TEST_EXE ${TEST_ID})
  set(HIBRIDON_MODULES_DIR ${hibridon_BINARY_DIR}) # location where to find the module files (*.mod)

  add_executable(${EXE_NAME} ${POT_SRC_FILE} ${hibridon_SOURCE_DIR}/src/himain.F90)
  target_compile_definitions(${EXE_NAME} PUBLIC T_MATRIX_SIZE=${p_T_MATRIX_SIZE})
  target_link_libraries(${EXE_NAME} hib ${HIBRIDON_LINK_LIBRARIES})
  target_include_directories(${EXE_NAME} PUBLIC ${hibridon_SOURCE_DIR}/src)
  target_include_directories(${EXE_NAME} PUBLIC ${HIBRIDON_MODULES_DIR})
  target_compile_options(${EXE_NAME} PRIVATE ${HIBRIDON_COMPILE_OPTIONS})

  if(ENABLE_CODE_COVERAGE)
    target_link_options(${EXE_NAME} PUBLIC ${GCC_COVERAGE_LINK_FLAGS})
  endif()
endfunction(add_hibexe)

file (GLOB userscripts ${CMAKE_CURRENT_BINARY_DIR}/*.user.cmake)
foreach(userscript IN LISTS userscripts)
   include (${userscript})
   add_hibexe("${EXE_NAME}" "${POT_SRC_FILE}" "${p_T_MATRIX_SIZE}")
   message ("Added ${userscript} user executable")
   unset(EXE_NAME)
   unset(POT_SRC_FILE)
   unset(p_T_MATRIX_SIZE)
endforeach()

