#------------------------------------------------------------------------------
# boh2 test
#------------------------------------------------------------------------------

# hib_boh2
add_executable(hib_boh2 src/pot/pot_boh2.F src/himain.F)
target_compile_definitions(hib_boh2 PUBLIC KMAX=120)
target_compile_definitions(hib_boh2 PUBLIC T_MATRIX_SIZE=kmax)
target_link_libraries(hib_boh2 hib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(hib_boh2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(hib_boh2 PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_compile_options(hib_boh2 PRIVATE ${HIBRIDON_COMPILE_OPTIONS})


add_test(NAME hib_boh2_test_setup
    COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_boh2; \
                    cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/boh2_test.com \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Boh2_bound.inp \
                        ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_boh2/"
)

add_test(NAME hib_boh2_build
    COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hib_boh2)

# $basedir/bin/progs/hib_boh2_151 <arn2_test.com
add_test(NAME hib_boh2_test_run
    COMMAND bash -c "cat ./boh2_test.com | ${CMAKE_CURRENT_BINARY_DIR}/hib_boh2"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests/hib_boh2"
)

# Check that outputs from test runs are the same as the ref ones
add_test(NAME hib_boh2_test_check
    COMMAND bash -c "./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Boh2_bou.evl ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_boh2/Boh2_bou.evl"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_test(NAME hib_boh2_test_cleanup
    COMMAND ls -R ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_boh2;
)


set_tests_properties(hib_boh2_test_setup PROPERTIES FIXTURES_SETUP hib_boh2_resources)
set_tests_properties(hib_boh2_test_cleanup PROPERTIES FIXTURES_CLEANUP hib_boh2_resources)

set_tests_properties(hib_boh2_test_run PROPERTIES FIXTURES_REQUIRED hib_boh2_resources)
set_tests_properties(hib_boh2_test_run PROPERTIES DEPENDS hib_boh2_build)
set_tests_properties(hib_boh2_test_check PROPERTIES DEPENDS check_outputs_build)
