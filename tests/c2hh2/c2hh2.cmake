#------------------------------------------------------------------------------
# C2H - oH2
#------------------------------------------------------------------------------

set(TEST_ID c2hh2)
set(TEST_POT_SRC_FILE "pot_c2hh2_12_6.F90")
set(TEST_POT_DATA_FILES "pot_c2hh2_12_6.dat")
set(TEST_COMMAND_FILE "test_c2hh2.com")
set(TEST_INPUT_FILES "C2hoh2.inp")
set(TEST_OUTPUT_FILES "Job1.xxsc")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")


#------------------------------------------------------------------------------
# C2H - oH2 (quick version for test coverage)
#------------------------------------------------------------------------------

set(TEST_ID c2hh2_quick)
set(TEST_POT_SRC_FILE "pot_c2hh2_12_6.F90")
set(TEST_POT_DATA_FILES "pot_c2hh2_12_6.dat")
set(TEST_COMMAND_FILE "test_c2hh2_quick.com")
set(TEST_INPUT_FILES "C2hoh2q.inp")
set(TEST_OUTPUT_FILES "C2hh2q1.xxsc")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")