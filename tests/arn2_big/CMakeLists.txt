#------------------------------------------------------------------------------
# arn2_big test
#------------------------------------------------------------------------------

set(TEST_ID arn2_big)
set(TEST_POT_SRC_FILE "pot_arn2.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "arn2_big.com")
set(TEST_INPUT_FILES "Arn2_test.inp")
set(TEST_OUTPUT_FILES "Csbtest1.ics Ccbrstest1.ics")
set(TEST_KMAX 36)
set(TEST_T_MATRIX_SIZE kbig)
set(TEST_LABELS "coverage")
add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
