#------------------------------------------------------------------------------
# h2ohe test
#------------------------------------------------------------------------------

# hib_h2ohe
add_executable(hib_h2ohe src/pot/pot_h2ohe.F src/himain.F)
target_compile_definitions(hib_h2ohe PUBLIC KMAX=45)
target_compile_definitions(hib_h2ohe PUBLIC T_MATRIX_SIZE=kmax)
target_link_libraries(hib_h2ohe hib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(hib_h2ohe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(hib_h2ohe PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_compile_options(hib_h2ohe PRIVATE ${HIBRIDON_COMPILE_OPTIONS})


add_test(NAME hib_h2ohe_test_setup
    COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_h2ohe; \
                    mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_h2ohe/potdata; \
                    cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/h2ohe_test.com \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/H2ohe_test.inp \
                        ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_h2ohe/ ; \
                    cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/h2o_coefd.dat \
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests/h2o_params.dat \
                    ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_h2ohe/potdata"
)

add_test(NAME hib_h2ohe_build
    COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hib_h2ohe)

# $basedir/bin/progs/hib_h2ohe_151 < h2ohe_test.com
add_test(NAME hib_h2ohe_test_run
    COMMAND bash -c "cat ./h2ohe_test.com | ${CMAKE_CURRENT_BINARY_DIR}/hib_h2ohe"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests/hib_h2ohe"
)

# Check that outputs from test runs are the same as the ref ones
add_test(NAME hib_h2ohe_test_check
    COMMAND bash -c "./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/H2ohe1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_h2ohe/H2ohe1.ics"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_test(NAME hib_h2ohe_test_cleanup
    COMMAND rm -R ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_h2ohe;
)


set_tests_properties(hib_h2ohe_test_setup PROPERTIES FIXTURES_SETUP hib_h2ohe_resources)
set_tests_properties(hib_h2ohe_test_cleanup PROPERTIES FIXTURES_CLEANUP hib_h2ohe_resources)

set_tests_properties(hib_h2ohe_test_run PROPERTIES FIXTURES_REQUIRED hib_h2ohe_resources)
set_tests_properties(hib_h2ohe_test_run PROPERTIES DEPENDS hib_h2ohe_build)
set_tests_properties(hib_h2ohe_test_check PROPERTIES DEPENDS check_outputs_build)
