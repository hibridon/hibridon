################################################################################
#
#    LIBRARY SOURCE FILES
#
################################################################################

set(ENABLE_ASSERTS "ON" CACHE BOOL "to enable or disable assert checks for debug builds (mostly for developers).")

set(LINK_TIME_OPTIMIZATION "ON" CACHE BOOL "to enable or disable link time optimizations on release builds for compilers that support it (slightly improves performance but greatly slows down link time)")

add_library(hib
    ancou.F90
    assert.f90
    grovec.F90
    hiamp.F90
    hibasis.F90
    hibasutil.F90
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
    hiprsbr.F90
    histmix.F90
    hisystem.F90
    hitensor.F90
    hitrnprt.F90
    hiutil.F90
    hivector.F90
    hipotutil.F90
    hiversion.F90
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

################################################################################
#
#    COMPILE OPTIONS
#
################################################################################

# GNU (gfortran)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(hib
      PUBLIC
      # Non-specific options
      -std=legacy                               # Allow pre-Fortran 77 syntax (e.g. arbitrary length arrays)
      -ffree-line-length-none                   # Allow arbitrary long lines. Needed as preprocessing could generate long line lengths.
      -fno-automatic                            # Put everything on the heap
#      -fdefault-integer-8
#      -finteger-4-integer-8
      $<$<PLATFORM_ID:Linux>:-march=native>     # Tune generated code for host's micro-architecture
      $<$<PLATFORM_ID:Darwin>:-mtune=native>    # Tune generated code for host's micro-architecture (march not working on Darwin)
      $<$<PLATFORM_ID:Linux>:-mcmodel=large>    # Required on Linux
      # Config-specific options: RELEASE
      $<$<CONFIG:RELEASE>:-O3>                  # Optimization level at 3 for Release
      $<$<AND:$<BOOL:${LINK_TIME_OPTIMIZATION}>,$<CONFIG:RELEASE>>:-flto>  # Activate link time optimizations so that gfortran can inline some function calls
      #$<$<CONFIG:RELEASE>:-fopt-info>           # You can get some information from gfortran with the flag -fopt-info that will tell you about optimizations missed or performed by the compiler
      $<$<CONFIG:RELEASE>:-finit-local-zero>    # Init variables to zero/false/null
      # Config-specific options: DEBUG
      $<$<CONFIG:DEBUG>:-O0>                    # Optimization level at 0
      $<$<CONFIG:DEBUG>:-g>                     # Include symbols in executable for easier debugging
      $<$<CONFIG:DEBUG>:-fno-omit-frame-pointer>
      $<$<CONFIG:DEBUG>:-fbacktrace>            # Generates extra information in the object file to provide source file traceback information when a severe error occurs at run time
      $<$<CONFIG:DEBUG>:-Wall>                  # Enable all warnings
      $<$<CONFIG:DEBUG>:-Wextra>                # Enable extra warnings
      $<$<CONFIG:DEBUG>:-fsanitize=address>     # Address sanitizer (beware: tests showed that with -fsanitize=address memory is not actually released when deallocate is called. Is this a true momory leak in address sanitizer?)
      $<$<CONFIG:DEBUG>:-Wuninitialized>        # Emit warnings for uninitialized variables 
      $<$<CONFIG:DEBUG>:-finit-integer=333333333> # initialize variables to non zero values to ease detection of uninitialized variables (some compilers simply don't initialize variables)
      $<$<CONFIG:DEBUG>:-finit-real=snan> # initialize variables to non zero values to ease detection of uninitialized variables (some compilers simply don't initialize variables)
      $<$<BOOL:${ENABLE_CODE_COVERAGE}>:--coverage> # Code coverage (same as -fprofile-arcs -ftest-coverage at compile time)
      $<$<BOOL:${ENABLE_PROFILING}>:-g>         # The profiler requires both the debug and profile directives (-g and -p)
      $<$<BOOL:${ENABLE_PROFILING}>:-p>         # The profiler requires both the debug and profile directives (-g and -p)
    )
# Intel (ifort)
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
    target_compile_options(hib
      PUBLIC
      # Non-specific options
      # force arrays to be allocated on the heap instead of the stack, this removes segmentation faults crash, not sure why
      # the following option means : automatic arrays and arrays created for temporary computations are allocated on the stack if their size can be determined at compile time and if it doesn't exceed 10kb.
      # this option seems necessary for big values of kmax (eg kmax=5000), according to bdesrousseaux, otherwise the user will experience a segmentation fault. I guess that without this option, the stack becomes too small to contain such big arrays...
      "-heap-arrays 10"                           # Put everything bigger than 10 bytes on the heap
      # "-extend-source 132"                        # Allow arbitrary long lines (132 seems the longest allowed with ifort). Needed as preprocessing could generate long lines.
      -no-wrap-margin                             # Don't wrap output files
      $<$<PLATFORM_ID:Linux>:-mcmodel=large>      # Required on Linux
      # Config-specific options: RELEASE
      $<$<CONFIG:RELEASE>:-O3>                    # Optimization level at 3 for Release
      $<$<CONFIG:RELEASE>:-xHost>                 # Generate instructions for the highest instruction set available on the compilation host processor
      $<$<AND:$<BOOL:${LINK_TIME_OPTIMIZATION}>,$<CONFIG:RELEASE>>:-ipo>  # activate interprocediral optimization (aka link time optimization)
      $<$<CONFIG:RELEASE>:-init=zero>             # Init variables to zero/false/null
      # Config-specific options: DEBUG
      $<$<CONFIG:DEBUG>:-O0>                      # Disable all optimizations
      $<$<CONFIG:DEBUG>:-g>                   # Generates complete debugging information
      $<$<CONFIG:DEBUG>:-traceback>               # Generates extra information in the object file to provide source file traceback information when a severe error occurs at run time
      $<$<CONFIG:DEBUG>:-fp-stack-check>          # Tell the compiler to generate extra code after every function call to ensure that the floating-point stack is in the expected state
      $<$<CONFIG:DEBUG>:-init=zero>             # Init integer and logical variables to zero/false/null instead of random to avoid random bugs
      $<$<CONFIG:DEBUG>:-init=snan>             # Init real variables to signaling nans to detect uninitialized variables
      $<$<CONFIG:DEBUG>:-init=array>            # also initialize arrays to avoid random behaviours caused by use of uninitialized variables
      #$<$<CONFIG:DEBUG>:-warn all>                # Enable all warnings
      $<$<BOOL:${ENABLE_PROFILING}>:-g>         # The profiler requires both the debug and profile directives (-g and -p)
      $<$<BOOL:${ENABLE_PROFILING}>:-p>         # The profiler requires both the debug and profile directives (-g and -p)
    )
endif()

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

# The following is for preprocessing of hiversion.F90
target_compile_definitions(hib PUBLIC _BUILD_VERS_=${PROJECT_VERSION})
target_compile_definitions(hib PUBLIC _BUILD_TYPE_="${CMAKE_BUILD_TYPE}")
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





