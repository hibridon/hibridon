# Function that configure executables for each test
function(add_hibridon_test TEST_ID TEST_POT_SRC_FILE TEST_POT_DATA_FILES TEST_COMMAND_FILE TEST_INPUT_FILES TEST_OUTPUT_FILES TEST_T_MATRIX_SIZE TEST_LABELS)

  set(TEST_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
  set(TEST_SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
  set(TEST_EXE ${TEST_ID})


  # Create an executable for that test
  add_hibexe(${TEST_EXE} ${TEST_SRC_DIR}/${TEST_POT_SRC_FILE} ${TEST_T_MATRIX_SIZE})

  file(MAKE_DIRECTORY ${TEST_BUILD_DIR}/potdata)
  separate_arguments(TEST_COMMAND_FILE)
  foreach(COMMAND_FILE ${TEST_COMMAND_FILE})
    configure_file(${TEST_SRC_DIR}/${COMMAND_FILE} ${TEST_BUILD_DIR}/${COMMAND_FILE} COPYONLY)
  endforeach()
  separate_arguments(TEST_INPUT_FILES)
  foreach(INPUT_FILE ${TEST_INPUT_FILES})
    configure_file(${TEST_SRC_DIR}/${INPUT_FILE} ${TEST_BUILD_DIR}/${INPUT_FILE} COPYONLY)
  endforeach()
  separate_arguments(TEST_POT_DATA_FILES)
  foreach(POT_DATA_FILE ${TEST_POT_DATA_FILES})
    configure_file(${TEST_SRC_DIR}/${POT_DATA_FILE} ${TEST_BUILD_DIR}/potdata/${POT_DATA_FILE} COPYONLY)
  endforeach()


  # Create a test that runs the executable using the provided com files
  add_test(NAME ${TEST_ID}_test_run COMMAND ${TEST_EXE} --kmax ${TEST_KMAX} --com ${TEST_COMMAND_FILE} )
  set_property(TEST ${TEST_ID}_test_run PROPERTY LABELS ${TEST_ID} ${TEST_LABELS})

  # Check the outputs of the test
  separate_arguments(TEST_OUTPUT_FILES)
  foreach(OUTPUT_FILE ${TEST_OUTPUT_FILES})
    add_test(NAME ${TEST_ID}_test_check_${OUTPUT_FILE} COMMAND check_outputs ${TEST_SRC_DIR}/${OUTPUT_FILE} ${TEST_BUILD_DIR}/${OUTPUT_FILE})
    set_property(TEST ${TEST_ID}_test_check_${OUTPUT_FILE} PROPERTY LABELS ${TEST_ID} ${TEST_LABELS})
  endforeach()


endfunction()