#------------------------------------------------------------------------------
# Ar - CH4
#------------------------------------------------------------------------------

set(TEST_ID arch4)
set(TEST_POT_SRC_FILE "pot_arch4.F")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "test_arch4.com")
set(TEST_INPUT_FILES "Arch4_a.inp")
set(TEST_OUTPUT_FILES "Job1.xxsc")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")