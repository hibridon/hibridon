#------------------------------------------------------------------------------
# Ar - CH4
#------------------------------------------------------------------------------

set(TEST_ID arch4)
set(TEST_POT_SRC_FILE "pot_arch4.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "test_arch4.com")
set(TEST_INPUT_FILES "Arch4_a.inp")
set(TEST_OUTPUT_FILES "Job1.xxsc")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# Ar - CH4 (quick version for test coverage)
#------------------------------------------------------------------------------

set(TEST_ID arch4_quick)
set(TEST_POT_SRC_FILE "pot_arch4.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "test_arch4_quick.com")
set(TEST_INPUT_FILES "Arch4_aq.inp")
set(TEST_OUTPUT_FILES "Arch4q1.xxsc Saveq1.inp")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage quick")
add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")