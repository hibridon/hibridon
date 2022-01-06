# Function that configure executables for each test
function(add_hibridon_test TEST_ID TEST_POT_SRC_FILE TEST_POT_DATA_FILES TEST_COMMAND_FILE TEST_INPUT_FILES TEST_OUTPUT_FILES TEST_T_MATRIX_SIZE)
  message("Adding test ${TEST_ID}")

  set(TEST_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
  set(TEST_SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
  set(TEST_EXE ${TEST_ID})


  # Create an executable for that test
  add_hibexe(${TEST_EXE} ${TEST_SRC_DIR}/${TEST_POT_SRC_FILE} ${TEST_T_MATRIX_SIZE})



  # Create a setup test that will copy input, com and potdata files to the binary dir
  # This should be replaced with add_custom_commands...
  add_test(NAME ${TEST_ID}_test_setup
  COMMAND bash -c "mkdir -p ${TEST_BUILD_DIR}/potdata && \
                  for input_file in ${TEST_INPUT_FILES} ${TEST_COMMAND_FILE}; \
                  do \
                      cp ${TEST_SRC_DIR}/$input_file ${TEST_BUILD_DIR}/ ; \
                  done ; \
                  for pot_data_file in ${TEST_POT_DATA_FILES}; \
                  do \
                      cp ${TEST_SRC_DIR}/$pot_data_file ${TEST_BUILD_DIR}/potdata ; \
                  done"
)
set_tests_properties(${TEST_ID}_test_setup PROPERTIES LABELS ${TEST_ID})

# Create a test that builds the test executable
#add_test(NAME ${TEST_ID}_build
#COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${TEST_EXE})

# Create a test that runs the executable using the provided com files
add_test(NAME ${TEST_ID}_test_run COMMAND ${TEST_EXE} --kmax ${TEST_KMAX} --com ${TEST_COMMAND_FILE} )
set_tests_properties(${TEST_ID}_test_run PROPERTIES LABELS ${TEST_ID})

# Check the outputs of the test
separate_arguments(TEST_OUTPUT_FILES)
foreach(OUTPUT_FILE ${TEST_OUTPUT_FILES})
  add_test(NAME ${TEST_ID}_test_check_${OUTPUT_FILE} COMMAND check_outputs ${TEST_SRC_DIR}/${OUTPUT_FILE} ${TEST_BUILD_DIR}/${OUTPUT_FILE})
  set_tests_properties(${TEST_ID}_test_check_${OUTPUT_FILE} PROPERTIES LABELS ${TEST_ID})
endforeach()

#add_test(NAME ${TEST_ID}_test_cleanup COMMAND ls -R ${TEST_BUILD_DIR};)

endfunction()
