#------------------------------------------------------------------------------
# arohx_ump4 test
#------------------------------------------------------------------------------

# hib_arohx_ump4
add_executable(hib_arohx_ump4 src/pot/pot_arohx_ump4.F src/himain.F)
target_compile_definitions(hib_arohx_ump4 PUBLIC KMAX=301)
target_compile_definitions(hib_arohx_ump4 PUBLIC T_MATRIX_SIZE=kmax)
target_link_libraries(hib_arohx_ump4 hib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(hib_arohx_ump4 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(hib_arohx_ump4 PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_compile_options(hib_arohx_ump4 PRIVATE ${HIBRIDON_COMPILE_OPTIONS})


add_test(NAME hib_arohx_ump4_test_setup
    COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arohx_ump4; \
                    cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/aroh_jatest.com \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Aroh_jatest.inp \
                        ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arohx_ump4/"
)

add_test(NAME hib_arohx_ump4_build
    COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hib_arohx_ump4)

# $basedir/bin/progs/hib_arohx_ump4_301 < h2ohe_test.com
add_test(NAME hib_arohx_ump4_test_run
    COMMAND bash -c "cat ./aroh_jatest.com | ${CMAKE_CURRENT_BINARY_DIR}/hib_arohx_ump4"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arohx_ump4"
)

# Check that outputs from test runs are the same as the ref ones
add_test(NAME hib_arohx_ump4_test_check
    COMMAND bash -c "./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Aroh_new1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arohx_ump4/Aroh_new1.ics"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_test(NAME hib_arohx_ump4_test_cleanup
    COMMAND ls -R ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_arohx_ump4;
)


set_tests_properties(hib_arohx_ump4_test_setup PROPERTIES FIXTURES_SETUP hib_arohx_ump4_resources)
set_tests_properties(hib_arohx_ump4_test_cleanup PROPERTIES FIXTURES_CLEANUP hib_arohx_ump4_resources)

set_tests_properties(hib_arohx_ump4_test_run PROPERTIES FIXTURES_REQUIRED hib_arohx_ump4_resources)
set_tests_properties(hib_arohx_ump4_test_run PROPERTIES DEPENDS hib_arohx_ump4_build)
set_tests_properties(hib_arohx_ump4_test_check PROPERTIES DEPENDS check_outputs_build)
