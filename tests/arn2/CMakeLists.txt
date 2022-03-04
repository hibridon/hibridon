#------------------------------------------------------------------------------
# arn2 test
#------------------------------------------------------------------------------

set(TEST_ID arn2)
set(TEST_POT_SRC_FILE "pot_arn2.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "arn2_test.com")
set(TEST_INPUT_FILES "Arn2_test.inp Arn2_dxsec.inp Arn2.fluxinp")
set(TEST_OUTPUT_FILES "Cctest1.ics Cstest1.ics Ccrstest1.ics")
set(TEST_KMAX 151)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")
add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
