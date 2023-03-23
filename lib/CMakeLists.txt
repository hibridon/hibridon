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
#      -fdefault-integer-8
#      -finteger-4-integer-8
      $<$<PLATFORM_ID:Linux>:-mcmodel=large>    # Required on Linux
      $<$<BOOL:${OpenMP_Fortran_FOUND}>:-fopenmp> # Compile with fopenmp option if OpenMP libs are found
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
      $<$<CONFIG:DEBUG>:-fsanitize=address>     # Address sanitizer
      $<$<AND:$<NOT:$<BOOL:${ENABLE_UNINIT_VAR_RUNTIME_DETECTOR}>>,$<CONFIG:DEBUG>>:-Wuninitialized>        # Emit warnings for uninitialized variables. Disable -Wuninitialized when ENABLE_UNINIT_VAR_RUNTIME_DETECTOR is on because the the -finit-* options then used make -Wuninitialized have no effect, see gfortran documentation :
      # Finally, note that enabling any of the -finit-* options will silence warnings that would have been emitted by -Wuninitialized for the affected local variables.

      # handle ENABLE_UNINIT_VAR_RUNTIME_DETECTOR
      # don't initialize local variables to zero : initialize them with dummy values to force failures which help detect uninitialized variables as -Wuninitialized doesn't detect everything (see issue #38 or https://stackoverflow.com/questions/39591300/gfortran-wuninitialized-flag-misses-variable-in-do-loop)
      # by default, gfortran initializes integers to 0, but not ifort : as a result, some bugs in the code are hidden with gfortran default options
      # set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fno-init-local-zero")
      # initialize variables to something else than 0 to force the programe to behave badly in case of unitialized variables
      $<$<AND:$<BOOL:${ENABLE_UNINIT_VAR_RUNTIME_DETECTOR}>,$<CONFIG:DEBUG>>:-finit-integer=333333333>
      $<$<AND:$<BOOL:${ENABLE_UNINIT_VAR_RUNTIME_DETECTOR}>,$<CONFIG:DEBUG>>:-finit-real=snan>
      $<$<AND:$<BOOL:${ENABLE_UNINIT_VAR_RUNTIME_DETECTOR}>,$<CONFIG:DEBUG>>:-ffpe-trap=invalid,zero,overflow>


      $<$<BOOL:${ENABLE_CODE_COVERAGE}>:--coverage -Wno-coverage-invalid-line-number> # Code coverage (same as -fprofile-arcs -ftest-coverage at compile time)
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
      -heap-arrays                                # Put everything on the heap
      # "-extend-source 132"                        # Allow arbitrary long lines (132 seems the longest allowed with ifort). Needed as preprocessing could generate long lines.
      -no-wrap-margin                             # Don't wrap output files
      $<$<PLATFORM_ID:Linux>:-mcmodel=large>      # Required on Linux
      $<$<BOOL:${OpenMP_Fortran_FOUND}>:-qopenmp> # Compile with qopenmp option if OpenMP libs are found
      # Config-specific options: RELEASE
      $<$<CONFIG:RELEASE>:-O3>                    # Optimization level at 3 for Release
      $<$<AND:$<BOOL:${LINK_TIME_OPTIMIZATION}>,$<CONFIG:RELEASE>>:-ipo>  # activate interprocediral optimization (aka link time optimization)
      $<$<CONFIG:RELEASE>:-init=zero>             # Init variables to zero/false/null
      # Config-specific options: DEBUG
      $<$<CONFIG:DEBUG>:-O0>                      # Disable all optimizations
      $<$<CONFIG:DEBUG>:-g>                   # Generates complete debugging information
      $<$<CONFIG:DEBUG>:-traceback>               # Generates extra information in the object file to provide source file traceback information when a severe error occurs at run time
      $<$<CONFIG:DEBUG>:-fp-stack-check>          # Tell the compiler to generate extra code after every function call to ensure that the floating-point stack is in the expected state
      #$<$<CONFIG:DEBUG>:-warn all>                # Enable all warnings

      # handle ENABLE_UNINIT_VAR_RUNTIME_DETECTOR
      # Force à NaN toutes les variables de type intrinsèque ainsi que les tableaux non initialisés. Cette option implique -fpe0. Pour éviter de récupérer des exceptions, qui ne soient pas relatives à des variables non-initialisées, nous recommandons de réduire le niveau d'optimisation à -O1 ou -O0 ou alors d'utiliser -fp-speculation=safe pour faire la détection de variables non-initialisés.
      $<$<AND:$<BOOL:${ENABLE_UNINIT_VAR_RUNTIME_DETECTOR}>,$<CONFIG:DEBUG>>:-init=zero>   # Init integer and logical variables to zero/false/null instead of random to avoid random bugs
      $<$<AND:$<BOOL:${ENABLE_UNINIT_VAR_RUNTIME_DETECTOR}>,$<CONFIG:DEBUG>>:-init=snan>    # Init real variables to signaling nans to detect uninitialized variables
      $<$<AND:$<BOOL:${ENABLE_UNINIT_VAR_RUNTIME_DETECTOR}>,$<CONFIG:DEBUG>>:-init=arrays>    # also initialize arrays to avoid random behaviours caused by use of uninitialized variables

      # check uses of uninitialized variables in run time
      # this is very useful as these are always bugs that cause the program to behave randomly with ifort (gfortran initializes data with zero)
      $<$<AND:$<BOOL:${ENABLE_UNINIT_VAR_RUNTIME_DETECTOR}>,$<CONFIG:DEBUG>>:-check uninit>

      # Cette combinaison d'options stoppe l'exécution dès qu'une exception (overflow, underflow, division par zéro, opération invalide,…) se produit; elle indique à quel niveau du code elle s'est produite. Le contrôle s'opère dans chaque subroutine, contrairement à l'option -fpe0 qui agit uniquement sur le programme principal. 
      # $<$<AND:$<BOOL:${ENABLE_UNINIT_VAR_RUNTIME_DETECTOR}>,$<CONFIG:DEBUG>>:-fpe-all=0>
      # $<$<AND:$<BOOL:${ENABLE_UNINIT_VAR_RUNTIME_DETECTOR}>,$<CONFIG:DEBUG>>:-no-ftz=0>
      # $<$<AND:$<BOOL:${ENABLE_UNINIT_VAR_RUNTIME_DETECTOR}>,$<CONFIG:DEBUG>>:-traceback=0>

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





