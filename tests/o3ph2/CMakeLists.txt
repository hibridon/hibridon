#------------------------------------------------------------------------------
# o3ph2 test
#------------------------------------------------------------------------------

set(TEST_ID o3ph2)
set(TEST_POT_SRC_FILE "pot_o3ph2.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "o3ph2.com")
set(TEST_INPUT_FILES "O3ph2.inp")
set(TEST_OUTPUT_FILES "Job1.ics Job1.xsc Job1.pcs Job1.dcs")
set(TEST_KMAX 30)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# o3ph2 test
#------------------------------------------------------------------------------

set(TEST_ID o3ph2_quick)
set(TEST_POT_SRC_FILE "pot_o3ph2.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "o3ph2_quick.com")
set(TEST_INPUT_FILES "O3ph2q.inp")
set(TEST_OUTPUT_FILES "Jobq1.ics Jobq1.xsc Jobq1.pcs Jobq1.dcs")
set(TEST_KMAX 30)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
