#------------------------------------------------------------------------------
# ohh2_j13 test
# OH--H2 j2=1,3 built-in basis
#------------------------------------------------------------------------------

set(TEST_ID ohh2_j13)
set(TEST_POT_SRC_FILE "pot_ohh2.F90")
set(TEST_POT_DATA_FILES "pot_ohh2_ccsdf12_avtzbf.dat")
set(TEST_COMMAND_FILE "test_ohh2.com")
set(TEST_INPUT_FILES "Ohh2.inp")
set(TEST_OUTPUT_FILES "Job1.ics Job1.pcs Job1.xxsc")
set(TEST_KMAX 1000)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# ohh2_j13 test (quick version for test coverage)
# OH--H2 j2=1,3 built-in basis
#------------------------------------------------------------------------------

set(TEST_ID ohh2_j13_quick)
set(TEST_POT_SRC_FILE "pot_ohh2.F90")
set(TEST_POT_DATA_FILES "pot_ohh2_ccsdf12_avtzbf.dat")
set(TEST_COMMAND_FILE "test_ohh2_quick.com")
set(TEST_INPUT_FILES "Ohh2q.inp")
set(TEST_OUTPUT_FILES "Ohh2q1.ics Ohh2q1.pcs Ohh2q1.xxsc")
set(TEST_KMAX 1000)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
