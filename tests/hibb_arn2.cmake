#------------------------------------------------------------------------------
#hibb_arn2 test
#------------------------------------------------------------------------------

add_executable(hibb_arn2 src/pot/pot_arn2.F src/himain.F)
target_compile_definitions(hibb_arn2 PUBLIC KMAX=36)
target_compile_definitions(hibb_arn2 PUBLIC T_MATRIX_SIZE=kbig)
target_link_libraries(hibb_arn2 hib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(hibb_arn2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(hibb_arn2 PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_compile_options(hibb_arn2 PRIVATE ${HIBRIDON_COMPILE_OPTIONS})


add_test(NAME hibb_arn2_test_setup
    COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/hibb_arn2; \
                    cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/arn2_big.com \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Arn2_test.inp \
                        ${CMAKE_CURRENT_BINARY_DIR}/tests/hibb_arn2/"
)

add_test(NAME hibb_arn2_build
    COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hibb_arn2)

# $basedir/bin/progs/hibb_arn2_36 <arn2_big.com
add_test(NAME hibb_arn2_test_run
    COMMAND bash -c "cat ./arn2_big.com | ${CMAKE_CURRENT_BINARY_DIR}/hibb_arn2"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests/hibb_arn2"
)

# Check that outputs from test runs are the same as the ref ones
add_test(NAME hibb_arn2_test_check
    COMMAND bash -c "./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Ccbtest1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/hibb_arn2/Ccbtest1.ics && \
    ./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Ccbrstest1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/hibb_arn2/Ccbrstest1.ics && \
    ./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Csbtest1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/hibb_arn2/Csbtest1.ics"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_test(NAME hibb_arn2_test_cleanup
    COMMAND rm -R ${CMAKE_CURRENT_BINARY_DIR}/tests/hibb_arn2;
)

set_tests_properties(hibb_arn2_test_setup PROPERTIES FIXTURES_SETUP hibb_arn2_resources)
set_tests_properties(hibb_arn2_test_cleanup PROPERTIES FIXTURES_CLEANUP hibb_arn2_resources)

set_tests_properties(hibb_arn2_test_run PROPERTIES FIXTURES_REQUIRED hibb_arn2_resources)
set_tests_properties(hibb_arn2_test_run PROPERTIES DEPENDS hibb_arn2_build)
set_tests_properties(hibb_arn2_test_check PROPERTIES DEPENDS check_outputs_build)
