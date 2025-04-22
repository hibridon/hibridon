################################################################################
#
#    LIBRARY SOURCE FILES
#
################################################################################

set(ENABLE_ASSERTS "ON" CACHE BOOL "to enable or disable assert checks for debug builds (mostly for developers).")

set(LINK_TIME_OPTIMIZATION "OFF" CACHE BOOL "to enable or disable link time optimizations on release builds for compilers that support it (slightly improves performance but greatly slows down link time)")

# With gfortran, variables are initialized to 0 by default, but it's not the case with ifort compiler. Therefore, all variables need to be initialized to guarantee the program works in a deterministic way on every compiler.
set(ENABLE_UNINIT_VAR_RUNTIME_DETECTOR "ON" CACHE BOOL "if set, causes the use of uninitilazed variables to trigger an exception on debug builds, thus making the program stop. This is useful to detect uninitialized memory issues.")

add_library(hib
    ancou.F90
    assert.f90
    grovec.F90
    hiamp.F90
    hibasis.F90
    hibasutil.F90
    hiblas.F90
    hibound.F90
    hibrid1.F90
    hibrid2.F90
    hibrid3.F90
    hibrid4.F90
    hibrid5.F90
    hicommon.F90
    hidifcrs.F90
    hiflow.F90
    hihypxsc.F90
    hiiolib_f.F90
    himatrix.F90
    himem.F90
    hinput.F90
    hiparcst.F90
    hipot.F90
    hiprsbr.F90
    hismat.F90
    histmix.F90
    hisystem.F90
    hitensor.F90
    hitrnprt.F90
    hitestptn.F90
    hiutil.F90
    hivector.F90
    hipotutil.F90
    hiversion.F90
    hiwav.F90
    parbas.F90
    parpot.F90
    bases/hiba01_1sg.F90
    bases/hiba02_2sg.F90
    bases/hiba03_2pi.F90
    bases/hiba04_sgpi.F90
    bases/hiba05_pi.F90
    bases/hiba06_stp.F90
    bases/hiba07_13p.F90
    bases/hiba08_2mol.F90
    bases/hiba09_stpln.F90
    bases/hiba10_22p.F90
    bases/hiba11_1del.F90
    bases/hiba12_h2p.F90
    bases/hiba13_h3p.F90
    bases/hiba14_2del.F90
    bases/hiba15_diat2p.F90
    bases/hiba16_astp.F90
    bases/hiba17_ch2x.F90
    bases/hiba18_stp1.F90
    bases/hiba19_sgpi1.F90
    bases/hiba20_2pi1sg.F90
    bases/hiba21_stp1sg.F90
    bases/hiba22_1d3p.F90
    bases/hiba23_3p2s.F90
    bases/hiba24_sphtp.F90
    bases/hiba25_1sg1sg.F90
    bases/hiba26_2sg1sg.F90
    bases/hiba27_astp1.F90
    bases/hiba28_3sg1sg.F90
    bases/hiba29_astp2.F90
    bases/hiba30_astp3.F90
    )

# This for preprocessed includes from the subdirectory bases to work
target_include_directories(hib PRIVATE ./)

# Link Blas and Lapack libraries
target_link_libraries(hib INTERFACE "${BLAS_LIBRARIES}")
target_link_libraries(hib INTERFACE "${LAPACK_LIBRARIES}")
# Link OpenMP if found
if(OpenMP_Fortran_FOUND)
    target_link_libraries(hib INTERFACE OpenMP::OpenMP_Fortran)
endif()

set_compile_options(hib "OFF")



################################################################################
#
#    LINK OPTIONS
#
################################################################################

# GNU (gfortran)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(hib
      PUBLIC
      $<$<CONFIG:DEBUG>:-fsanitize=address>         # Address sanitizer (beware: tests showed that with -fsanitize=address memory is not actually released when deallocate is called. Is this a true momory leak in address sanitizer?)
      $<$<BOOL:${ENABLE_CODE_COVERAGE}>:--coverage> # Code coverage (same as -lgcov at link time)
      $<$<BOOL:${ENABLE_PROFILING}>:-g>         # The profiler requires both the debug and profile directives (-g and -p)
      $<$<BOOL:${ENABLE_PROFILING}>:-p>         # The profiler requires both the debug and profile directives (-g and -p)
      $<$<AND:$<BOOL:${LINK_TIME_OPTIMIZATION}>,$<CONFIG:RELEASE>>:-flto>  # Activate link time optimizations so that gfortran can inline some function calls
      #$<$<CONFIG:RELEASE>:-fopt-info>           # You can get some information from gfortran with the flag -fopt-info that will tell you about optimizations missed or performed by the compiler
    )
# Intel (ifort)
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
    target_link_libraries(hib
      PUBLIC
      $<$<BOOL:${ENABLE_PROFILING}>:-g>         # The profiler requires both the debug and profile directives (-g and -p)
      $<$<BOOL:${ENABLE_PROFILING}>:-p>         # The profiler requires both the debug and profile directives (-g and -p)
      $<$<AND:$<BOOL:${LINK_TIME_OPTIMIZATION}>,$<CONFIG:RELEASE>>:-ipo>  # activate interprocediral optimization (aka link time optimization)
    )
endif()


################################################################################
#
#    COMPILE DEFINITIONS
#
################################################################################

# The following is to define the machine type and compiler
target_compile_definitions(hib PUBLIC HIB_UNIX_X86=1)
target_compile_definitions(hib PUBLIC HIB_UNIX=1)
if (CMake_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_definitions(hib PUBLIC HIB_UNIX_GFORTRAN=1)
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
  target_compile_definitions(hib PUBLIC HIB_UNIX_IFORT=1)
endif()

if(ENABLE_CODE_COVERAGE)
  if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  else()
    message(FATAL_ERROR "ENABLE_CODE_COVERAGE is not yet supported with compiler ${Fortran_COMPILER_NAME}")
  endif()

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


  # build 
  add_custom_target(clean_coverages VERBATIM COMMAND bash -c "\
    for info_file in $(find ${CMAKE_BINARY_DIR} -name \"*.info\") ; \
    do \
      echo deleting $info_file ; \
      rm $info_file ; \
    done ;"
    COMMENT "deleting the *.info *.gcda files found in ${CMAKE_BINARY_DIR} to erase all traces of previous exe runs for code coverages")

  add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/coverage/total.info COMMAND bash -c "\
    mkdir -p ${CMAKE_BINARY_DIR}/coverage ; \
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


# make the code aware that ENABLE_UNINIT_VAR_RUNTIME_DETECTOR is on, to allow the code to use workarounds to prevent the code from stopping because of use on unintialized variables in some libraries we have no control (eg mkl's dsyevr function uses uninitialized variables) 
target_compile_definitions(hib PUBLIC $<$<AND:$<BOOL:${ENABLE_UNINIT_VAR_RUNTIME_DETECTOR}>,$<CONFIG:DEBUG>>:FLOATING_POINT_ERRORS_CAUSE_EXCEPTIONS=1>)

# The following is for preprocessing of hiversion.F90
target_compile_definitions(hib PUBLIC _BUILD_VERS_=${PROJECT_VERSION})
target_compile_definitions(hib PUBLIC _BUILD_TYPE_="${CMAKE_BUILD_TYPE}")
execute_process (COMMAND bash -c "git -C ${CMAKE_CURRENT_LIST_DIR} rev-parse HEAD 2>/dev/null" OUTPUT_VARIABLE GIT_REVISION)
if(NOT GIT_REVISION STREQUAL "")
  string(STRIP ${GIT_REVISION} GIT_REVISION)
endif()
target_compile_definitions(hib PUBLIC _GIT_REVISION_="${GIT_REVISION}")
target_compile_definitions(hib PUBLIC _BUILD_SYS1_="${CMAKE_HOST_SYSTEM_NAME}")
target_compile_definitions(hib PUBLIC _BUILD_SYS2_="${CMAKE_HOST_SYSTEM_VERSION}")
target_compile_definitions(hib PUBLIC _BUILD_PROC_="${CMAKE_HOST_SYSTEM_PROCESSOR}")
target_compile_definitions(hib PUBLIC _BUILD_COMP_="${CMAKE_Fortran_COMPILER}")

# The following is for preprocessing of hiutil.F90
target_compile_definitions(hib PUBLIC _HELPDIR_="${CMAKE_CURRENT_SOURCE_DIR}/../doc/")
target_compile_definitions(hib PUBLIC _BUILD_FLAGS_="")#"${HIBRIDON_COMPILE_OPTIONS}")
string(TIMESTAMP TODAY "%Y-%m-%d at %H:%M")
target_compile_definitions(hib PUBLIC _BUILD_DATE_="${TODAY}")
# The following is for preprocessing of hiutil.F90
target_compile_definitions(hib PUBLIC _HELPDIR_="${CMAKE_CURRENT_SOURCE_DIR}/../doc/")
# The following is to limit the max jtot number
target_compile_definitions(hib PUBLIC MAX_NJTOT=1000)
target_compile_definitions(hib PUBLIC "$<$<CONFIG:Release>:DISABLE_HIB_ASSERT=1>")
# asserts can be activated on non-release builds (eg debug build)
if(NOT ENABLE_ASSERTS)
  target_compile_definitions(hib PUBLIC "$<$<CONFIG:Debug>:DISABLE_HIB_ASSERT=1>")
endif()





