#------------------------------------------------------------------------------
# ca3phe test
#------------------------------------------------------------------------------

set(TEST_ID ca3phe)
set(TEST_POT_SRC_FILE "pot_ca3phe.F")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "ca3phe_test.com")
set(TEST_INPUT_FILES "Ca3phe.inp")
set(TEST_OUTPUT_FILES "")
set(TEST_KMAX 81)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")
