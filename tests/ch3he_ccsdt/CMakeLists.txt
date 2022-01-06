#------------------------------------------------------------------------------
# ch3he_ccsdt test
#------------------------------------------------------------------------------

set(TEST_ID ch3he_ccsdt)
set(TEST_POT_SRC_FILE "pot_ch3he_ccsdt.F90")
set(TEST_POT_DATA_FILES "ch3he_pot.dat")
set(TEST_COMMAND_FILE "ch3he_test.com")
set(TEST_INPUT_FILES "Ch3he_test.inp")
set(TEST_OUTPUT_FILES "Ch3he1.ics")
set(TEST_KMAX 151)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# ch3he_ccsdt test (Quick verison)
#------------------------------------------------------------------------------

set(TEST_ID ch3he_ccsdt_quick)
set(TEST_POT_SRC_FILE "pot_ch3he_ccsdt.F90")
set(TEST_POT_DATA_FILES "ch3he_pot.dat")
set(TEST_COMMAND_FILE "ch3he_test_quick.com")
set(TEST_INPUT_FILES "Ch3he_testq.inp")
set(TEST_OUTPUT_FILES "Ch3heq1.ics")
set(TEST_KMAX 151)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
