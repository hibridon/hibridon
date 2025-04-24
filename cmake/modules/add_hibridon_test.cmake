# Function that configure executables for each test
function(add_hibridon_test TEST_ID TEST_POT_SRC_FILE TEST_POT_DATA_FILES TEST_COMMAND_FILE TEST_INPUT_FILES TEST_OUTPUT_FILES TEST_T_MATRIX_SIZE TEST_LABELS)
  # TEST_ID: eg "arn2"
  # TEST_POT_SRC_FILE: eg "pot_c2hh2_12_6.F90"
  # TEST_POT_DATA_FILES: eg "pot_c2hh2_12_6.dat"
  # TEST_COMMAND_FILE: eg "test_c2hh2.com"
  # TEST_INPUT_FILES: eg "Arn2_test.inp Arn2_dxsec.inp Arn2.fluxinp"
  # TEST_OUTPUT_FILES: eg "Cctest1.ics Cstest1.ics Ccrstest1.ics"
  # TEST_T_MATRIX_SIZE: eg 151
  # TEST_LABELS: eg "coverage quick"
  separate_arguments(TEST_LABELS)
  set(TEST_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
  set(TEST_SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
  get_filename_component(TEST_EXE ${TEST_POT_SRC_FILE} NAME_WE)
  
  # Create an executable for that test
  if (NOT TARGET ${TEST_EXE})
    add_hibexe(${TEST_EXE} ${TEST_SRC_DIR}/${TEST_POT_SRC_FILE} ${TEST_T_MATRIX_SIZE})
  endif()

  # create a cmake test that creates the directories for the hibridon test
  add_test(NAME ${TEST_ID}_test_setup
    COMMAND bash -c "mkdir -p ${TEST_BUILD_DIR}/potdata && \
                    mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/coverage && \
                    for input_file in ${TEST_INPUT_FILES} ${TEST_COMMAND_FILE}; \
                    do \
                        cp ${TEST_SRC_DIR}/$input_file ${TEST_BUILD_DIR}/ ; \
                    done ; \
                    for pot_data_file in ${TEST_POT_DATA_FILES}; \
                    do \
                        cp ${TEST_SRC_DIR}/$pot_data_file ${TEST_BUILD_DIR}/potdata ; \
                    done"
  )
  set_property(TEST ${TEST_ID}_test_setup PROPERTY LABELS ${TEST_ID} ${TEST_LABELS})
  
  separate_arguments(TEST_OUTPUT_FILES)

  if(TEST_OUTPUT_FILES)  # only add this if the list of output files is not empty
    # remove output files before running the test as hibridon appends to some existing files, and this would cause the output file checker to fail the test
    add_test(NAME ${TEST_ID}_test_prolog COMMAND ${CMAKE_COMMAND} -E remove ${TEST_OUTPUT_FILES})
    set_property(TEST ${TEST_ID}_test_prolog PROPERTY WORKING_DIRECTORY ${TEST_BUILD_DIR})
    set_property(TEST ${TEST_ID}_test_prolog PROPERTY LABELS ${TEST_ID} ${TEST_LABELS})
  endif()

  # Create a cmake test that builds the hibridon test
  add_test(NAME ${TEST_ID}_build
  COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${TEST_EXE})
  set_property(TEST ${TEST_ID}_build PROPERTY LABELS ${TEST_ID} ${TEST_LABELS})

  # Create a test that runs the executable using the provided com files
  add_test(NAME ${TEST_ID}_test_run COMMAND ${TEST_EXE} --kmax ${TEST_KMAX} --com ${TEST_COMMAND_FILE} )
  set_property(TEST ${TEST_ID}_test_run PROPERTY LABELS ${TEST_ID} ${TEST_LABELS})
  set_property(TEST ${TEST_ID}_test_run PROPERTY DEPENDS ${TEST_ID}_test_prolog)

  # Check the outputs of the test
  foreach(OUTPUT_FILE ${TEST_OUTPUT_FILES})
    add_test(NAME ${TEST_ID}_test_check_${OUTPUT_FILE} COMMAND check_outputs ${TEST_SRC_DIR}/${OUTPUT_FILE} ${TEST_BUILD_DIR}/${OUTPUT_FILE})
    set_property(TEST ${TEST_ID}_test_check_${OUTPUT_FILE} PROPERTY LABELS ${TEST_ID} ${TEST_LABELS})
  endforeach()

  if(ENABLE_CODE_COVERAGE)
    add_test(NAME ${TEST_ID}_save_coverage
      COMMAND lcov -c -d ${CMAKE_BINARY_DIR}/ -o ${CMAKE_CURRENT_BINARY_DIR}/coverage/${TEST_ID}.info;
    )
    set_property(TEST ${TEST_ID}_save_coverage PROPERTY LABELS ${TEST_ID} ${TEST_LABELS})

    add_test(NAME ${TEST_ID}_cleanup_coverage
      COMMAND find ${CMAKE_BINARY_DIR} -name "*.gcda" -delete;
    )
    set_property(TEST ${TEST_ID}_cleanup_coverage PROPERTY LABELS ${TEST_ID} ${TEST_LABELS})
  endif()

  if(GENERATE_PROFILING_PDF)
    add_test(NAME ${TEST_ID}_build_profiling_pdf
      COMMAND  bash -c "${PROFILING_GPROF_EXE} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_EXE} > ${TEST_ID}_gprofout.txt && cat ${TEST_ID}_gprofout.txt | ${PROFILING_GPROF2DOT_EXE} | ${PROFILING_DOT_EXE} -Tpdf -o ${TEST_ID}_call_graph.pdf"
      WORKING_DIRECTORY "${TEST_BUILD_DIR}")
    set_property(TEST ${TEST_ID}_build_profiling_pdf PROPERTY LABELS ${TEST_ID} ${TEST_LABELS})
  endif()

  add_test(NAME ${TEST_ID}_test_cleanup
    COMMAND ls -R ${TEST_BUILD_DIR};
  )
  set_property(TEST ${TEST_ID}_test_cleanup PROPERTY LABELS ${TEST_ID} ${TEST_LABELS})
  
  # test dependencies
  set_tests_properties(${TEST_ID}_test_setup PROPERTIES FIXTURES_SETUP ${TEST_ID}_resources)
  if(ENABLE_CODE_COVERAGE)
    set_tests_properties(${TEST_ID}_cleanup_coverage PROPERTIES FIXTURES_CLEANUP ${TEST_ID}_resources)
  endif()
  set_tests_properties(${TEST_ID}_test_cleanup PROPERTIES FIXTURES_CLEANUP ${TEST_ID}_resources)

  set_tests_properties(${TEST_ID}_test_run PROPERTIES FIXTURES_REQUIRED ${TEST_ID}_resources)
  set_tests_properties(${TEST_ID}_test_run PROPERTIES DEPENDS ${TEST_ID}_build)
  # set_tests_properties(${TEST_ID}_test_check PROPERTIES DEPENDS check_outputs_build)
  foreach(OUTPUT_FILE ${TEST_OUTPUT_FILES})
    set_tests_properties("${TEST_ID}_test_check_${OUTPUT_FILE}" PROPERTIES DEPENDS check_outputs_build)
  endforeach()

  if(ENABLE_CODE_COVERAGE)
    # set_tests_properties(${TEST_ID}_save_coverage PROPERTIES DEPENDS ${TEST_ID}_test_check)
    foreach(OUTPUT_FILE ${TEST_OUTPUT_FILES})
      set_tests_properties(${TEST_ID}_save_coverage PROPERTIES DEPENDS "${TEST_ID}_test_check_${OUTPUT_FILE}")
    endforeach()
    set_tests_properties(${TEST_ID}_cleanup_coverage PROPERTIES DEPENDS ${TEST_ID}_save_coverage)
  endif()
 
  if(GENERATE_PROFILING_PDF)
    set_tests_properties(${TEST_ID}_build_profiling_pdf PROPERTIES LABELS ${TEST_ID} ${TEST_LABELS})
  endif()


endfunction()
