#------------------------------------------------------------------------------
# arn2 test
#------------------------------------------------------------------------------

# hib_arn2
add_executable(hib_arn2 src/pot/pot_arn2.F src/himain.F)
target_compile_definitions(hib_arn2 PUBLIC KMAX=151)
target_compile_definitions(hib_arn2 PUBLIC T_MATRIX_SIZE=kmax)
target_link_libraries(hib_arn2 hib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(hib_arn2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(hib_arn2 PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_compile_options(hib_arn2 PRIVATE ${HIBRIDON_COMPILE_OPTIONS})


add_test(NAME hib_arn2_test_setup
    COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arn2; \
                    cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/arn2_test.com \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Arn2_test.inp \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Arn2_dxsec.inp \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Arn2.fluxinp \
                        ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arn2/"
)

add_test(NAME hib_arn2_build
    COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hib_arn2)

# $basedir/bin/progs/hib_arn2_151 <arn2_test.com
# todo: fix memory leak and remove export ASAN_OPTIONS=detect_leaks=0
add_test(NAME hib_arn2_test_run
    COMMAND bash -c "export ASAN_OPTIONS=detect_leaks=0; cat ./arn2_test.com | ${CMAKE_CURRENT_BINARY_DIR}/hib_arn2"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arn2"
)

# Check that outputs from test runs are the same as the ref ones
add_test(NAME hib_arn2_test_check
    COMMAND bash -c "./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Cctest1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arn2/Cctest1.ics && \
    ./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Ccrstest1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arn2/Ccrstest1.ics && \
    ./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Cstest1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arn2/Cstest1.ics"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_test(NAME hib_arn2_test_cleanup
    COMMAND rm -R ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arn2;
)


set_tests_properties(hib_arn2_test_setup PROPERTIES FIXTURES_SETUP hib_arn2_resources)
set_tests_properties(hib_arn2_test_cleanup PROPERTIES FIXTURES_CLEANUP hib_arn2_resources)

set_tests_properties(hib_arn2_test_run PROPERTIES FIXTURES_REQUIRED hib_arn2_resources)
set_tests_properties(hib_arn2_test_run PROPERTIES DEPENDS hib_arn2_build)
set_tests_properties(hib_arn2_test_check PROPERTIES DEPENDS check_outputs_build)
