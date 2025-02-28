# The following is to build hibridon with an user-supplied potential
function(add_hibexe EXE_NAME POT_SRC_FILE p_T_MATRIX_SIZE)

  message("Adding executable ${EXE_NAME} using potential ${POT_SRC_FILE}")
  set(HIBRIDON_MODULES_DIR ${hibridon_BINARY_DIR}) # location where to find the module files (*.mod)

  add_executable(${EXE_NAME} ${POT_SRC_FILE} ${hibridon_SOURCE_DIR}/src/himain.F90)
  
  target_link_libraries(${EXE_NAME} PUBLIC hib)
  
  target_compile_definitions(${EXE_NAME} PRIVATE T_MATRIX_SIZE=${p_T_MATRIX_SIZE})

  target_include_directories(${EXE_NAME} PRIVATE ${hibridon_SOURCE_DIR}/lib) # to find the included common files
  target_include_directories(${EXE_NAME} PRIVATE ${hibridon_BINARY_DIR}/lib) # to find the .mod files


################################################################################
#
#    COMPILE OPTIONS
#
################################################################################

# GNU (gfortran)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${EXE_NAME}
      PUBLIC
      # Non-specific options
      -std=legacy                               # Allow pre-Fortran 77 syntax (e.g. arbitrary length arrays)
      -fno-automatic                            # Put everything on the heap
      $<$<PLATFORM_ID:Linux>:-mcmodel=large>    # Required on Linux
      # Config-specific options: RELEASE
      $<$<CONFIG:RELEASE>:-O3>                  # Optimization level at 3 for Release
      $<$<CONFIG:RELEASE>:-finit-local-zero>    # Init variables to zero/false/null
      # Config-specific options: DEBUG
      $<$<CONFIG:DEBUG>:-Og>                    # Optimization level at 0 and generates complete debugging information
      $<$<CONFIG:DEBUG>:-fbacktrace>            # Generates extra information in the object file to provide source file traceback information when a severe error occurs at run time
      $<$<CONFIG:DEBUG>:-Wall>                  # Enable all warnings
      $<$<CONFIG:DEBUG>:-Wextra>                # Enable extra warnings
      $<$<CONFIG:DEBUG>:-fsanitize=address>     # Address sanitizer
      $<$<CONFIG:DEBUG>:-Wuninitialized>        # Emit warnings for uninitialized variables 
      $<$<BOOL:${ENABLE_CODE_COVERAGE}>:--coverage> # Code coverage (same as -fprofile-arcs -ftest-coverage at compile time)
    )
# Intel (ifort)
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
    target_compile_options(${EXE_NAME}
      PUBLIC
      # Non-specific options
      -heap-arrays                                # Put everything on the heap
      -no-wrap-margin                             # Don't wrap output files
      # Config-specific options: RELEASE
      $<$<CONFIG:RELEASE>:-O3>                    # Optimization level at 3 for Release
      $<$<CONFIG:RELEASE>:-init=zero>             # Init variables to zero/false/null
      # Config-specific options: DEBUG
      $<$<CONFIG:DEBUG>:-O0>                      # Disable all optimizations
      $<$<CONFIG:DEBUG>:-g>                       # Generates complete debugging information
      $<$<CONFIG:DEBUG>:-traceback>               # Generates extra information in the object file to provide source file traceback information when a severe error occurs at run time
      $<$<CONFIG:DEBUG>:-warn all>                # Enable all warnings
    )
endif()

################################################################################
#
#    LINK OPTIONS
#
################################################################################

# GNU (gfortran)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_link_options(${EXE_NAME}
      PUBLIC
      $<$<CONFIG:DEBUG>:-fsanitize=address>         # Address sanitizer
      $<$<BOOL:${ENABLE_CODE_COVERAGE}>:--coverage> # Code coverage (same as -lgcov at link time)
    )
# Intel (ifort)
endif()


# Install instructions
install(TARGETS ${EXE_NAME} DESTINATION bin)

endfunction(add_hibexe)
