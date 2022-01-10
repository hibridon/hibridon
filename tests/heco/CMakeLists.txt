#------------------------------------------------------------------------------
# heco test
#------------------------------------------------------------------------------

set(TEST_ID heco)
set(TEST_POT_SRC_FILE "pot_heco_sapt.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "heco_test.com")
set(TEST_INPUT_FILES "Heco_test.inp")
set(TEST_OUTPUT_FILES "Heco1.ics")
set(TEST_KMAX 55)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

