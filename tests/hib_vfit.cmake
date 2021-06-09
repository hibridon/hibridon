#------------------------------------------------------------------------------
# vfit test
#------------------------------------------------------------------------------

# hib_vfit
add_executable(hib_vfit src/pot/pot_vfit.F src/himain.F)
target_compile_definitions(hib_vfit PUBLIC KMAX=36)
target_compile_definitions(hib_vfit PUBLIC T_MATRIX_SIZE=kmax)
target_link_libraries(hib_vfit hib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(hib_vfit PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(hib_vfit PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_compile_options(hib_vfit PRIVATE ${HIBRIDON_COMPILE_OPTIONS})


add_test(NAME hib_vfit_test_setup
    COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_vfit/potdata \
                    &&  cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/vfit_test.com \
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/N2phetest.inp \
                        ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_vfit/ \
                    && cp ${CMAKE_CURRENT_SOURCE_DIR}/tests/FOLLMEG.BIN \
                        ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_vfit/potdata"
)

add_test(NAME hib_vfit_build
    COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hib_vfit)

# $basedir/bin/progs/hib_vfit_151 <vfit_test.com
add_test(NAME hib_vfit_test_run
    COMMAND bash -c "cat ./vfit_test.com | ${CMAKE_CURRENT_BINARY_DIR}/hib_vfit"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests/hib_vfit"
)

# Check that outputs from test runs are the same as the ref ones
add_test(NAME hib_vfit_test_check
    COMMAND bash -c "./check_outputs ${CMAKE_CURRENT_SOURCE_DIR}/tests/Vfit_test1.ics ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_vfit/Vfit_test1.ics"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_test(NAME hib_vfit_test_cleanup
    COMMAND ls -R ${CMAKE_CURRENT_BINARY_DIR}/tests/hib_vfit;
)

set_tests_properties(hib_vfit_test_setup PROPERTIES FIXTURES_SETUP hib_vfit_resources)
set_tests_properties(hib_vfit_test_cleanup PROPERTIES FIXTURES_CLEANUP hib_vfit_resources)

set_tests_properties(hib_vfit_test_run PROPERTIES FIXTURES_REQUIRED hib_vfit_resources)
set_tests_properties(hib_vfit_test_run PROPERTIES DEPENDS hib_vfit_build)
set_tests_properties(hib_vfit_test_check PROPERTIES DEPENDS check_outputs_build)
