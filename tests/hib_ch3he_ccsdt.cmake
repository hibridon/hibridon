#------------------------------------------------------------------------------
# ch3he_ccsdt test
#------------------------------------------------------------------------------

# hib_ch3he_ccsdt
add_executable(hib_ch3he_ccsdt src/pot/pot_ch3he_ccsdt.F src/himain.F)
target_compile_definitions(hib_ch3he_ccsdt PUBLIC KMAX=151)
target_compile_definitions(hib_ch3he_ccsdt PUBLIC T_MATRIX_SIZE=kmax)
target_link_libraries(hib_ch3he_ccsdt hib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(hib_ch3he_ccsdt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(hib_ch3he_ccsdt PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_compile_options(hib_ch3he_ccsdt PRIVATE ${HIBRIDON_COMPILE_OPTIONS})


add_test(NAME hib_ch3he_ccsdt_test_setup
    COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3he_ccsdt; \
                    mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3he_ccsdt/potdata; \
                    cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/ch3he_test.com \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Ch3he_test.inp \
                        ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3he_ccsdt/ ; \
                    cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/ch3he_pot.dat \
                       ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3he_ccsdt/potdata"
)

add_test(NAME hib_ch3he_ccsdt_build
    COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hib_ch3he_ccsdt)

# $basedir/bin/progs/hib_ch3he_ccsdt_151 < ch3he_test.com
add_test(NAME hib_ch3he_ccsdt_test_run
    COMMAND bash -c "cat ./ch3he_test.com | ${CMAKE_CURRENT_BINARY_DIR}/hib_ch3he_ccsdt"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3he_ccsdt"
)

# Check that outputs from test runs are the same as the ref ones
add_test(NAME hib_ch3he_ccsdt_test_check
    COMMAND bash -c "./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Ch3he1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3he_ccsdt/Ch3he1.ics"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_test(NAME hib_ch3he_ccsdt_test_cleanup
    COMMAND rm -R ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3he_ccsdt;
)


set_tests_properties(hib_ch3he_ccsdt_test_setup PROPERTIES FIXTURES_SETUP hib_ch3he_ccsdt_resources)
set_tests_properties(hib_ch3he_ccsdt_test_cleanup PROPERTIES FIXTURES_CLEANUP hib_ch3he_ccsdt_resources)

set_tests_properties(hib_ch3he_ccsdt_test_run PROPERTIES FIXTURES_REQUIRED hib_ch3he_ccsdt_resources)
set_tests_properties(hib_ch3he_ccsdt_test_run PROPERTIES DEPENDS hib_ch3he_ccsdt_build)
set_tests_properties(hib_ch3he_ccsdt_test_check PROPERTIES DEPENDS check_outputs_build)
