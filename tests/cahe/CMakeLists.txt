#------------------------------------------------------------------------------
# cahe test
#------------------------------------------------------------------------------

set(TEST_ID cahe)
set(TEST_POT_SRC_FILE "pot_cahe.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "cahe_test.com")
set(TEST_INPUT_FILES "Cahe_test.inp")
set(TEST_OUTPUT_FILES "Job1.ics")
set(TEST_KMAX 191)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# nologdprop test (test case where start==rendld)
#------------------------------------------------------------------------------

set(TEST_ID nologdprop)
set(TEST_POT_SRC_FILE "pot_cahe.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "nologdprop_test.com")
set(TEST_INPUT_FILES "Nologdprop_test.inp")
set(TEST_OUTPUT_FILES "Nologdprop1.ics")
set(TEST_KMAX 191)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
