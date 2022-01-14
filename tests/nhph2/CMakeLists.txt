#------------------------------------------------------------------------------
# nhph2 test
#------------------------------------------------------------------------------

set(TEST_ID nhph2)
set(TEST_POT_SRC_FILE "pot_nhh2_10_4.F90")
set(TEST_POT_DATA_FILES "pot_nhh2_10_4.dat")
set(TEST_COMMAND_FILE "nhph2.com")
set(TEST_INPUT_FILES "Nhph2.inp")
set(TEST_OUTPUT_FILES "Job1.xsc")
set(TEST_KMAX 250)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

set(TEST_ID nhph2_quick)
set(TEST_POT_SRC_FILE "pot_nhh2_10_4.F90")
set(TEST_POT_DATA_FILES "pot_nhh2_10_4.dat")
set(TEST_COMMAND_FILE "nhph2_quick.com")
set(TEST_INPUT_FILES "Nhph2q.inp")
set(TEST_OUTPUT_FILES "Jobq1.xsc")
set(TEST_KMAX 220)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
