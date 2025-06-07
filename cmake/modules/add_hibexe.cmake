# The following is to build hibridon with an user-supplied potential
function(add_hibexe EXE_NAME POT_SRC_FILE p_T_MATRIX_SIZE)

  message("Adding executable ${EXE_NAME} using potential ${POT_SRC_FILE}")
  set(HIBRIDON_MODULES_DIR ${hibridon_BINARY_DIR}) # location where to find the module files (*.mod)

  add_executable(${EXE_NAME} ${POT_SRC_FILE} ${hibridon_SOURCE_DIR}/src/himain.F90)
  
  target_link_libraries(${EXE_NAME} PUBLIC hib)
  
  target_compile_definitions(${EXE_NAME} PRIVATE T_MATRIX_SIZE=${p_T_MATRIX_SIZE})

  target_include_directories(${EXE_NAME} PRIVATE ${hibridon_SOURCE_DIR}/lib) # to find the included common files
  target_include_directories(${EXE_NAME} PRIVATE ${hibridon_BINARY_DIR}/lib) # to find the .mod files

  set_compile_options(${EXE_NAME} "ON")


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
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM") # Intel (ifx)
    target_link_options(${EXE_NAME}
      PUBLIC
      $<$<AND:$<BOOL:${ENABLE_UNINIT_VAR_RUNTIME_DETECTOR}>,$<CONFIG:DEBUG>>:-check uninit>  # https://community.intel.com/t5/Intel-Fortran-Compiler/Linking-errors-when-using-memory-sanitizer-in-fortran-project/td-p/1521476: When you compile with -check uninit (or -check all) you also need to link with that compiler option.
    )
endif()


# Install instructions
install(TARGETS ${EXE_NAME} DESTINATION bin)

endfunction(add_hibexe)
