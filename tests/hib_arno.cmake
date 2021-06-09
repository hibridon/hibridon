#------------------------------------------------------------------------------
# arno test
#------------------------------------------------------------------------------

add_executable(hib_arno src/pot/pot_arno.F src/himain.F)
target_compile_definitions(hib_arno PUBLIC KMAX=191)
target_compile_definitions(hib_arno PUBLIC T_MATRIX_SIZE=kmax)
target_link_libraries(hib_arno hib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(hib_arno PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(hib_arno PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_compile_options(hib_arno PRIVATE ${HIBRIDON_COMPILE_OPTIONS})


add_test(NAME hib_arno_test_setup
    COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arno \
                    &&  cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/arno_test.com \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Arno_test.inp \
                        ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arno"
)

add_test(NAME hib_arno_build
    COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hib_arno)

# $basedir/bin/progs/hib_arno_151 <arno_test.com
add_test(NAME hib_arno_test_run
    COMMAND bash -c "cat ./arno_test.com | ${CMAKE_CURRENT_BINARY_DIR}/hib_arno"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arno"
)

# Check that outputs from test runs are the same as the ref ones
add_test(NAME hib_arno_test_check
    COMMAND bash -c "./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Arno_tes1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arno/Arno_tes1.ics"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_test(NAME hib_arno_test_cleanup
    COMMAND ls -R ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arno;
)

set_tests_properties(hib_arno_test_setup PROPERTIES FIXTURES_SETUP hib_arno_resources)
set_tests_properties(hib_arno_test_cleanup PROPERTIES FIXTURES_CLEANUP hib_arno_resources)

set_tests_properties(hib_arno_test_run PROPERTIES FIXTURES_REQUIRED hib_arno_resources)
set_tests_properties(hib_arno_test_run PROPERTIES DEPENDS hib_arno_build)
set_tests_properties(hib_arno_test_check PROPERTIES DEPENDS check_outputs_build)
