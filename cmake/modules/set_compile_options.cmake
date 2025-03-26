
function(set_compile_options TARGET SAVE)
# TARGET: the cmake target identifier (eg the name of a library target such as hib) of the target that compile options need to be set
# SAVE: if true, all local variables receive the save attribute. This argument should be removed once all programs work without this option (using local variables with an implicit save attribute should be avoided as it's error prone)
################################################################################
#
#    COMPILE OPTIONS
#
################################################################################
# GNU (gfortran)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${TARGET}
      PUBLIC
      # Non-specific options
      -std=legacy                               # Allow pre-Fortran 77 syntax (e.g. arbitrary length arrays)
      -ffree-line-length-none                   # Allow arbitrary long lines. Needed as preprocessing could generate long line lengths.
#      -fdefault-integer-8
#      -finteger-4-integer-8
      $<$<BOOL:${SAVE}>:-fno-automatic>          # Put local variables on the heap instead of the stack. As a side effect, all affected variables receive the save attribute (persistance between calls). This flag conflicts with -fopenmp (f951: Warning: Flag ‘-fno-automatic’ overwrites ‘-frecursive’ implied by ‘-fopenmp’)
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
      $<$<CONFIG:DEBUG>:-Wno-conversion>        # Disable warnings such as Warning: Possible change of value in conversion from REAL(8) to INTEGER(4) at (1) [-Wconversion]
      $<$<CONFIG:DEBUG>:-Wno-compare-reals>     # Disable warnings such as "Warning: Equality comparison for REAL(8) at (1) [-Wcompare-reals]" as they are often fine when comparing with 2 hardcoded zero constants
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
    target_compile_options(${TARGET}
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
      $<$<CONFIG:DEBUG>:-g>                       # Generates complete debugging information
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
      $<$<CONFIG:DEBUG>:-warn all>                # Enable all warnings
    )
endif()

endfunction()