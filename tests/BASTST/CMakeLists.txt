#------------------------------------------------------------------------------
# Dummy tests
#------------------------------------------------------------------------------


set(TEST_POT_SRC_FILE "pot_bastst.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")

# BASTST for basis 01
set(TEST_ID bastst_01)
set(TEST_INPUT_FILES  "bastst_01.inp")
set(TEST_OUTPUT_FILES "Bastst_01.out")
set(TEST_COMMAND_FILE "bastst_01.com")
add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

# BASTST for basis 02
set(TEST_ID bastst_02)
set(TEST_INPUT_FILES  "bastst_02.inp")
set(TEST_OUTPUT_FILES "Bastst_02.out")
set(TEST_COMMAND_FILE "bastst_02.com")
add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
