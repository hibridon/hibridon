# Function that configure executables for each test
function(add_hibridon_test TEST_ID TEST_POT_SRC_FILE TEST_POT_DATA_FILES TEST_COMMAND_FILE TEST_INPUT_FILES TEST_OUTPUT_FILES TEST_T_MATRIX_SIZE)
  message("Adding test ${TEST_ID}")


  # Create an executable for that test
  add_hibexe(${TEST_ID} ${CMAKE_CURRENT_LIST_DIR}/${TEST_POT_SRC_FILE} ${TEST_T_MATRIX_SIZE})
  
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


endfunction()