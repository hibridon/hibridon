# Function that configure executables for each test
function(add_hibridon_test TEST_ID TEST_POT_SRC_FILE TEST_POT_DATA_FILES TEST_COMMAND_FILE TEST_INPUT_FILES TEST_OUTPUT_FILES TEST_T_MATRIX_SIZE)
  message("Adding test ${TEST_ID}")


  # Create an executable for that test
  add_hibexe(${TEST_ID} ${CMAKE_CURRENT_LIST_DIR}/${TEST_POT_SRC_FILE} ${TEST_T_MATRIX_SIZE})
  

  set(TEST_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/tests/${TEST_ID})
  
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

  #add_test(NAME ${TEST_ID}_build COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${TEST_EXE})
  
  add_test(NAME ${TEST_ID}_test_run
    COMMAND bash -c "cat ./${TEST_COMMAND_FILE} | ${CMAKE_CURRENT_BINARY_DIR}/${TEST_ID} --kmax ${TEST_KMAX} | tee ./${TEST_ID}.stdout"
    WORKING_DIRECTORY "${TEST_BUILD_DIR}"
  )

  add_test(NAME ${TEST_ID}_test_check
      COMMAND bash -c "set -o errexit; \
                        for output_file in ${TEST_OUTPUT_FILES} ; \
                        do \
                          ${hibridon_BINARY_DIR}/check_outputs ${TEST_SRC_DIR}/$output_file ${TEST_BUILD_DIR}/$output_file ; 
                        done"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )

endfunction()
