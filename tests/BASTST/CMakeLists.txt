#------------------------------------------------------------------------------
# Dummy tests
#------------------------------------------------------------------------------

set(TEST_ID bastst)
set(TEST_POT_SRC_FILE "pot_bastst.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "bastst.com")
set(TEST_INPUT_FILES "bastst_01.inp bastst_02.inp")
set(TEST_OUTPUT_FILES "Bastst_01.out Bastst_02.out")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")
add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
