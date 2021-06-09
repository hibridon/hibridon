#------------------------------------------------------------------------------
# heco_sapt test
#------------------------------------------------------------------------------

# hib_heco_sapt
add_executable(hib_heco_sapt src/pot/pot_heco_sapt.F src/himain.F)
target_compile_definitions(hib_heco_sapt PUBLIC KMAX=55)
target_compile_definitions(hib_heco_sapt PUBLIC T_MATRIX_SIZE=kmax)
target_link_libraries(hib_heco_sapt hib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(hib_heco_sapt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(hib_heco_sapt PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_compile_options(hib_heco_sapt PRIVATE ${HIBRIDON_COMPILE_OPTIONS})


add_test(NAME hib_heco_sapt_test_setup
    COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_heco_sapt; \
                    cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/heco_test.com \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/Heco_test.inp \
                        ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_heco_sapt/"
)

add_test(NAME hib_heco_sapt_build
    COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hib_heco_sapt)

# $basedir/bin/progs/hib_heco_sapt_151 < h2ohe_test.com
add_test(NAME hib_heco_sapt_test_run
    COMMAND bash -c "cat ./heco_test.com | ${CMAKE_CURRENT_BINARY_DIR}/hib_heco_sapt"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests/hib_heco_sapt"
)

# Check that outputs from test runs are the same as the ref ones
add_test(NAME hib_heco_sapt_test_check
    COMMAND bash -c "./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Heco1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_heco_sapt/Heco1.ics"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_test(NAME hib_heco_sapt_test_cleanup
    COMMAND rm -R ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_heco_sapt;
)


set_tests_properties(hib_heco_sapt_test_setup PROPERTIES FIXTURES_SETUP hib_heco_sapt_resources)
set_tests_properties(hib_heco_sapt_test_cleanup PROPERTIES FIXTURES_CLEANUP hib_heco_sapt_resources)

set_tests_properties(hib_heco_sapt_test_run PROPERTIES FIXTURES_REQUIRED hib_heco_sapt_resources)
set_tests_properties(hib_heco_sapt_test_run PROPERTIES DEPENDS hib_heco_sapt_build)
set_tests_properties(hib_heco_sapt_test_check PROPERTIES DEPENDS check_outputs_build)
