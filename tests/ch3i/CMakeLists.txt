#------------------------------------------------------------------------------
# ch3i test
#------------------------------------------------------------------------------

set(TEST_ID ch3i)
set(TEST_POT_SRC_FILE "pot_ch3i.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "ch3i_test.com")
set(TEST_INPUT_FILES "Ch3i.inp")
set(TEST_OUTPUT_FILES "Ch3itest.flx")
set(TEST_KMAX 100)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")