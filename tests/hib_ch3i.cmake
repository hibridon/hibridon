#------------------------------------------------------------------------------
# ch3i test
#------------------------------------------------------------------------------

add_executable(hib_ch3i src/pot/pot_ch3i.F src/himain.F)
target_compile_definitions(hib_ch3i PUBLIC KMAX=100)
target_compile_definitions(hib_ch3i PUBLIC T_MATRIX_SIZE=kmax)
target_link_libraries(hib_ch3i hib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(hib_ch3i PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(hib_ch3i PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_compile_options(hib_ch3i PRIVATE ${HIBRIDON_COMPILE_OPTIONS})


add_test(NAME hib_ch3i_test_setup
    COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3i \
                    &&  cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/ch3i_test.com \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Ch3i.inp \
                        ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3i"
)

add_test(NAME hib_ch3i_build
    COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hib_ch3i)

add_test(NAME hib_ch3i_test_run
    COMMAND bash -c "cat ./ch3i_test.com | ${CMAKE_CURRENT_BINARY_DIR}/hib_ch3i"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3i"
)

# Check that outputs from test runs are the same as the ref ones
# todo: replace with a better comparer that uses a tolerance (at the moment, we compate the meaningful text but not the last line because the numbers are so tiny that they differ from a build and another)
add_test(NAME hib_ch3i_test_check
    COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3i/flx_compare \
        && cat ${CMAKE_CURRENT_SOURCE_DIR}/tests/Ch3itest.flx | tail -35 | head -34 > ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3i/flx_compare/gold.txt \
        && cat ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3i/Ch3itest.flx | tail -35 | head -34 > ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3i/flx_compare/measured.txt \
        && diff ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3i/flx_compare/gold.txt ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3i/flx_compare/measured.txt"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_test(NAME hib_ch3i_test_cleanup
    COMMAND ls -R ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_ch3i;
)

set_tests_properties(hib_ch3i_test_setup PROPERTIES FIXTURES_SETUP hib_ch3i_resources)
set_tests_properties(hib_ch3i_test_cleanup PROPERTIES FIXTURES_CLEANUP hib_ch3i_resources)

set_tests_properties(hib_ch3i_test_run PROPERTIES FIXTURES_REQUIRED hib_ch3i_resources)
set_tests_properties(hib_ch3i_test_run PROPERTIES DEPENDS hib_ch3i_build)
set_tests_properties(hib_ch3i_test_check PROPERTIES DEPENDS check_outputs_build)
