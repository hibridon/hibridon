#------------------------------------------------------------------------------
# hecn_dgels test
#------------------------------------------------------------------------------

add_executable(hib_hecn_dgels src/pot/pot_hecn_dgels.F src/himain.F)
target_compile_definitions(hib_hecn_dgels PUBLIC KMAX=81)
target_compile_definitions(hib_hecn_dgels PUBLIC T_MATRIX_SIZE=kmax)
target_link_libraries(hib_hecn_dgels hib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(hib_hecn_dgels PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(hib_hecn_dgels PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_compile_options(hib_hecn_dgels PRIVATE ${HIBRIDON_COMPILE_OPTIONS})

add_test(NAME hib_hecn_dgels_test_setup
    COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_hecn_dgels/potdata \
                    &&  cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/hecnx_hyp.com\
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Hecnx.inp \
                        ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_hecn_dgels \
                    && cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/hecn_fitmlv.dat  \
                    ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_hecn_dgels/potdata"
)

add_test(NAME hib_hecn_dgels_build
    COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hib_hecn_dgels)

add_test(NAME hib_hecn_dgels_test_run
    COMMAND bash -c "cat ./hecnx_hyp.com | ${CMAKE_CURRENT_BINARY_DIR}/hib_hecn_dgels"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests/hib_hecn_dgels"
)

# Check that outputs from test runs are the same as the ref ones
# Hecn1.hfx
# Hecn1.ics
# Hecn1.pcs
# Hecn1.smt
# Hecn1.xsc
# Hecn.sav
# hecnx_hyp.com
# Hecnx.inp
# Outpt
add_test(NAME hib_hecn_dgels_test_check
    COMMAND bash -c "./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Hecn1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_hecn_dgels/Hecn1.ics"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_test(NAME hib_hecn_dgels_test_cleanup
    COMMAND ls -R ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_hecn_dgels;
)

set_tests_properties(hib_hecn_dgels_test_setup PROPERTIES FIXTURES_SETUP hib_hecn_dgels_resources)
set_tests_properties(hib_hecn_dgels_test_cleanup PROPERTIES FIXTURES_CLEANUP hib_hecn_dgels_resources)

set_tests_properties(hib_hecn_dgels_test_run PROPERTIES FIXTURES_REQUIRED hib_hecn_dgels_resources)
set_tests_properties(hib_hecn_dgels_test_run PROPERTIES DEPENDS hib_hecn_dgels_build)
set_tests_properties(hib_hecn_dgels_test_check PROPERTIES DEPENDS check_outputs_build)
