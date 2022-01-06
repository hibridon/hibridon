#------------------------------------------------------------------------------
# ohh2_j1 test
# OH--H2 j2=1 built-in basis
#------------------------------------------------------------------------------

set(TEST_ID ohh2_j1)
set(TEST_POT_SRC_FILE "pot_ohh2.F90")
set(TEST_POT_DATA_FILES "pot_ohh2_ccsdf12_avtzbf.dat")
set(TEST_COMMAND_FILE "test_ohh2.com")
set(TEST_INPUT_FILES "Ohh2.inp")
set(TEST_OUTPUT_FILES "Job1.ics Job1.pcs Job1.xxsc")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# ohh2_j1 test (quick version for code coverage)
# OH--H2 j2=1 built-in basis
#------------------------------------------------------------------------------

set(TEST_ID ohh2_j1_quick)
set(TEST_POT_SRC_FILE "pot_ohh2.F90")
set(TEST_POT_DATA_FILES "pot_ohh2_ccsdf12_avtzbf.dat")
set(TEST_COMMAND_FILE "ohh2_j1_quick.com")
set(TEST_INPUT_FILES "Ohh2_quick.inp")
set(TEST_OUTPUT_FILES "Quick1.ics Quick1.pcs Quick1.xxsc")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
