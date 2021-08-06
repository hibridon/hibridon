#------------------------------------------------------------------------------
# ohh2_j1 test
# OH--H2 j2=1 built-in basis
#------------------------------------------------------------------------------

set(TEST_ID ohh2_j1)
set(TEST_POT_SRC_FILE "pot_ohh2.F")
set(TEST_POT_DATA_FILES "pot_ohh2_ccsdf12_avtzbf.dat")
set(TEST_COMMAND_FILE "test_ohh2.com")
set(TEST_INPUT_FILES "Ohh2.inp")
set(TEST_OUTPUT_FILES "Job1.ics Job1.pcs Job1.xxsc")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")
