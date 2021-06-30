#------------------------------------------------------------------------------
# boh2 test
#------------------------------------------------------------------------------

set(TEST_ID boh2)
set(TEST_POT_SRC_FILE "pot_boh2.F")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "boh2_test.com")
set(TEST_INPUT_FILES "Boh2_bound.inp")
set(TEST_OUTPUT_FILES "Boh2_bou.evl")
set(TEST_KMAX 120)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")

