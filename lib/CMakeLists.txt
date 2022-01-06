################################################################################
#
#    LIBRARY SOURCE FILES
#
################################################################################


add_library(hib
    assert.f90
    hiamp.F90
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


################################################################################
#
#    COMPILE OPTIONS
#
################################################################################

# GNU (gfortran)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(hib
    PUBLIC
    -std=legacy
    -fmax-stack-var-size=1
    $<$<AND:UNIX,$<NOT:APPLE>>: -mcmodel=large>
    $<$<CONFIG:RELEASE>:-finit-local-zero> # Initializes variables and arrays to 0
    $<$<CONFIG:RELEASE>:-O0>
    $<$<CONFIG:DEBUG>:-Og>
    $<$<CONFIG:DEBUG>:-Wall>
    $<$<CONFIG:DEBUG>:-Wextra>
    $<$<CONFIG:DEBUG>:-ffpe-trap=invalid,zero,overflow>
    )
# Intel (ifort)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    # TODO
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
elseif (CMake_Fortran_COMPILER_ID STREQUAL "Intel")
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

# The ffollowing is for preprocessing of hhelp.F90 and hiutil.F90
target_compile_definitions(hib PUBLIC _HELPDIR_="${CMAKE_CURRENT_SOURCE_DIR}/doc/")
target_compile_definitions(hib PUBLIC _BUILD_FLAGS_="")#"${HIBRIDON_COMPILE_OPTIONS}")
string(TIMESTAMP TODAY "%Y-%m-%d at %H:%M")
target_compile_definitions(hib PUBLIC _BUILD_DATE_="${TODAY}")
# The ffollowing is for preprocessing of hhelp.F90 and hiutil.F90
target_compile_definitions(hib PUBLIC _HELPDIR_="${CMAKE_CURRENT_SOURCE_DIR}/doc/")
# The following is to limit the max jtot number
target_compile_definitions(hib PUBLIC MAX_NJTOT=1000)
target_compile_definitions(hib PUBLIC DISABLE_HIB_ASSERT=1)


# Link Blas and Lapack libraries
target_link_libraries(hib INTERFACE "${BLAS_LIBRARIES}")
target_link_libraries(hib INTERFACE "${LAPACK_LIBRARIES}")


